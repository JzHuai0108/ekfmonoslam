#!/usr/bin/env python
# Expect time stamped transformations and gyro angular rates
# Output time offset between the two signal streams

from __future__ import print_function
import argparse
import bisect
import csv
import os
import subprocess
import sys


import numpy as np
from scipy import signal
from time_alignment_plotting_tools import (plot_results,
                                           plot_angular_velocities)

def read_time_stamped_angular_rate_from_file(csv_file):
  """
  Reads time stamped poses from a CSV or TXT file.
  Assumes the following line format:
    timestamp [s], gx [rad/s], gy [rad/s], gz [rad/s], ...
  """
  delimiter = ','
  with open(csv_file, 'r') as csvfile:
    lineswithcomma = 0
    totallines = 0
    for line in csvfile:
        if ',' in line:
            lineswithcomma += 1
        totallines += 1
    if lineswithcomma / float(totallines) < 0.1:
        delimiter = ' '

  with open(csv_file, 'r') as csvfile:
    csv_reader = csv.reader(csvfile, delimiter=delimiter, quotechar='|')
    datalist = list()
    for row in csv_reader:
        if '#' in row[0] or '%' in row[0] or '[' in row[0]:
            continue
        datalist.append(row)
    time_stamped_rates = np.array(datalist)
    time_stamped_rates = time_stamped_rates.astype(float)
    time_stamped_rates = time_stamped_rates[:, :4]

  times = time_stamped_rates[:, 0].copy()
  rates = time_stamped_rates[:, 1:]
  return (time_stamped_rates.copy(), times, rates)


class FilteringConfig:

  def __init__(self):
    self.smoothing_kernel_size_A = 25
    self.clipping_percentile_A = 99.5
    self.smoothing_kernel_size_B = 25
    self.clipping_percentile_B = 99.0


def filter_and_smooth_angular_velocity(angular_velocity,
                                       low_pass_kernel_size, clip_percentile, plot=False):
  """Reduce the noise in a velocity signal."""

  max_value = np.percentile(angular_velocity, clip_percentile)
  print("Clipping angular velocity norms to {} rad/s ...".format(max_value))
  angular_velocity_clipped = np.clip(angular_velocity, -max_value, max_value)
  print("Done clipping angular velocity norms...")

  low_pass_kernel = np.ones((low_pass_kernel_size, 1)) / low_pass_kernel_size
  print("Smoothing with kernel size {} samples...".format(low_pass_kernel_size))

  angular_velocity_smoothed = signal.correlate(angular_velocity_clipped,
                                               low_pass_kernel, 'same')

  print("Done smoothing angular velocity norms...")

  if plot:
    plot_angular_velocities("Angular Velocities", angular_velocity,
                            angular_velocity_smoothed, True)

  return angular_velocity_smoothed.copy()


def calculate_time_offset_from_signals(times_A, signal_A,
                                       times_B, signal_B,
                                       plot=False, block=True):
  """ Calculates the time offset between signal A and signal B. """
  convoluted_signals = signal.correlate(signal_B, signal_A)
  dt_A = np.mean(np.diff(times_A))
  offset_indices = np.arange(-len(signal_A) + 1, len(signal_B))
  max_index = np.argmax(convoluted_signals)
  offset_index = offset_indices[max_index]

  # quadratic interpolation
  valm1 = convoluted_signals[max_index - 1]
  valmax = convoluted_signals[max_index]
  valp1 = convoluted_signals[max_index + 1]
  coeffs = [(valp1 + valm1 - 2*valmax)/2, (valp1 - valm1)/2, valmax]

  offset_fitted = - coeffs[1] / (2*coeffs[0])
  valmax_fitted = - coeffs[1]*coeffs[1] / (4*coeffs[0]) + coeffs[2]
  if offset_fitted >= -1 and offset_fitted <= 1:
      print('Updated discrete lag index %d to %.3f and maxcorr %.4f to %.4f'
            ' by fit2' % (offset_index, offset_fitted + offset_index, valmax, valmax_fitted))
      offset_index += offset_fitted

  time_offset = dt_A * offset_index + times_B[0] - times_A[0]

  if plot:
    plot_results(times_A, times_B, signal_A, signal_B, convoluted_signals,
                 time_offset, block=block)
  return time_offset


def compute_angular_velocity_norms(angular_velocity, smoothing_kernel_size, clipping_percentile, plot=False):
  """angular velocity NX3"""
  angular_velocity_norms = []
  angular_velocity_size = len(angular_velocity)

  filter_and_smooth = True
  if filter_and_smooth:
    angular_velocity_filtered = filter_and_smooth_angular_velocity(
        angular_velocity, smoothing_kernel_size, clipping_percentile, plot)
  else:
    angular_velocity_filtered = angular_velocity.copy()

  for i in range(0, angular_velocity_size):
    angular_velocity_norms.append(
        np.linalg.norm(angular_velocity_filtered[i, :]))

  assert len(angular_velocity_norms) == angular_velocity_size
  return angular_velocity_norms

def resample_rates(times, rates, dt):
    """
    resample rates at the times specified at samples with cubit fit
    :return: resampled rates and samples
    """
    samples = np.arange(times[0], times[-1], dt)
    interpolation_needed = False
    for index, sample in enumerate(samples):
        if np.abs(sample - times[index]) > 1e-8:
            interpolation_needed = True
            break
    if not interpolation_needed:
        print("Interpolation is unneeded for times beginning at {} of median"
              " sample interval {}".format(times[0], np.median(np.diff(times))))
        return rates, times
    from scipy.interpolate import CubicSpline
    cs = CubicSpline(times, rates) # if rates is 3XN, then axis=1 has to be passed
    return cs(samples), samples

def calculate_time_offset(times_A, rates_A, times_B, rates_B, filtering_config, plot=False):
  """
  Calculate the time offset between the stamped rates_A and rates_B.
  take the norm and then apply a convolution to compute the best time alignment for the two sets of angular rates.
  Note that the accuracy of the time alignment is limited by the higher frequency
  of the two signals, i.e. by the smallest time interval between two rates.
  """

  # Get the two mean time steps. Take the smaller one for the interpolation.
  dt_A = np.median(np.diff(times_A))
  dt_B = np.median(np.diff(times_B))
  if dt_A >= dt_B:
    dt = dt_B
  else:
    dt = dt_A
  print("Interpolation interval {}".format(dt))

  rates_A_interp, samples_A = resample_rates(times_A, rates_A, dt)
  rates_B_interp, samples_B = resample_rates(times_B, rates_B, dt)

  # Compute angular velocity norms for the resampled rates.
  angular_velocity_norms_A = compute_angular_velocity_norms(
      rates_A_interp,
      filtering_config.smoothing_kernel_size_A,
      filtering_config.clipping_percentile_A, plot)
  angular_velocity_norms_B = compute_angular_velocity_norms(
      rates_B_interp,
      filtering_config.smoothing_kernel_size_B,
      filtering_config.clipping_percentile_B, plot)

  assert len(samples_A) == len(angular_velocity_norms_A)

  # Comput time offset.
  time_offset = calculate_time_offset_from_signals(
      samples_A, angular_velocity_norms_A, samples_B,
      angular_velocity_norms_B, plot, True)


  return time_offset, dt


def write_double_numpy_array_to_csv_file(array, csv_file):
  np.savetxt(csv_file, array, delimiter=", ", fmt="%.18f")

if __name__ == '__main__':
    """
    Perform time alignment between timestamped set of poses and set of angular rates
  
    The CSV file for poses, T_W_E, should have the following line format:
      timestamp [s], x [m], y [m], z [m], qx, qy, qz, qw
    If a TXT file is provided, the delimiter will be replaced by white spaces
    The CSV file for angular rate, \omega_{IS}^S, should have the following line format:
      timestamp [s], wx [rad/s], wy [rad/s], wz [rad/s], ...
    If a TXT file is provided, the delimiter will be replaced by white spaces
  
    """

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        '--angular_rate_file',
        required=False,
        help='The file containing the time stamped angular rates. (e.g. measured by gyro)')
    parser.add_argument(
        '--poses_B_H_file',
        required=False,
        help='The file containing the time stamped poses '
             '(e.g. hand poses in an earth-fixed Base frame).'
             ' You may either provide the angular_rate_file or poses_B_H_file.')
    parser.add_argument(
        '--poses_W_E_file',
        required=True,
        help='The file containing the time stamped poses. (e.g. Eye poses in an earth-fixed World frame)')

    parser.add_argument(
        '--upsample_se3_binary',
        required=True,
        help='The binary executable to upsample SE3 poses. usually under build folder generated by cmake and make')

    parser.add_argument(
        '--time_offset_output_csv_file', type=str,
        help='Write estimated time offset to this file in spatial-extrinsics csv format')
    parser.add_argument('--visualize',
                        help='Visualize the poses.',
                        action='store_true')
    if len(sys.argv) < 2:
        print("Example usage 1: {} --angular_rate_file "
              "/visual_inertial_data/8-4-2018/raw_accel_gyro.csv"
              " --poses_W_E_file /visual_inertial_data/8-4-2018/FrameTrajectory.txt"
              " --upsample_se3_binary /ekfmonoslam/SE3Interpolation/build/upsample_se3"
              " --visualize".format(sys.argv[0]))
        print("Example usage 2: {} --poses_B_H_file "
              "/visual_inertial_data/8-4-2018/lidarposes.txt"
              " --poses_W_E_file /visual_inertial_data/8-4-2018/FrameTrajectory.txt"
              " --upsample_se3_binary /ekfmonoslam/SE3Interpolation/build/upsample_se3"
              " --visualize".format(sys.argv[0]))
    args = parser.parse_args()

    print("Sampling angular rates from input SE3 poses of {}".format(args.poses_W_E_file))

    sampled_rate_csv = os.path.join(os.path.dirname(args.poses_W_E_file), "UpsampledPseudoImuE.csv")
    output_freq = 400
    commandlist = [args.upsample_se3_binary, args.poses_W_E_file, sampled_rate_csv, str(output_freq)]
    subprocess.call(commandlist)

    if args.poses_B_H_file:
        if args.angular_rate_file is not None:
            raise Exception("Only angular_rate_file or poses_B_H_file is needed")
        print("Sampling angular rates from input SE3 poses of {}".format(args.poses_B_H_file))
        args.angular_rate_file = os.path.join(os.path.dirname(args.poses_B_H_file), "UpsampledPseudoImuH.csv")
        output_freq = 400
        commandlist = [args.upsample_se3_binary, args.poses_B_H_file, args.angular_rate_file, str(output_freq)]
        subprocess.call(commandlist)

    print("Reading data from files...")
    (time_stamped_rates_E, times_E, rates_E) = read_time_stamped_angular_rate_from_file(sampled_rate_csv)
    print("Found ", time_stamped_rates_E.shape[
        0], " angular rates in file: ", sampled_rate_csv)

    (time_stamped_rates_H, times_H,
     rates_H) = read_time_stamped_angular_rate_from_file(args.angular_rate_file)
    print("Found ", time_stamped_rates_H.shape[0], " rates in file: ", args.angular_rate_file)

    print("Computing time offset...")
    filtering_config = FilteringConfig()
    filtering_config.visualize = args.visualize
    # TODO(mfehr): get filtering config from args!

    time_offset, dt_sigma = calculate_time_offset(times_H, rates_H, times_E,
                                        rates_E, filtering_config,
                                        filtering_config.visualize)

    print("Final time offset: ", time_offset, "+/-", dt_sigma, "s")
    print("time of rate_E - time_offset = time of rate_E in H clock")

    if args.time_offset_output_csv_file is not None:
        print("Writing time_offset to %s." % args.time_offset_output_csv_file)
        write_double_numpy_array_to_csv_file(np.array((time_offset,)),
                                             args.time_offset_output_csv_file)
