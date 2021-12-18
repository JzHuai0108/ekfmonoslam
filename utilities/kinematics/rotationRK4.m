% Huai 12/31/2014
% $\dot(q)_b^a=1/2A(\omega_{ab}^b)q_b^a$ is the differential equation,
% input omega contains timestamps and $\omega_{ab}^b$, q_k is $q_b^a$ at
% the first timestamp of omega, $q_b^a=[w;x;y;z]$.
function q_next_k = rotationRK4( omega, q_k)
times=omega(1,:);
dt=times(2:end)-times(1:end-1);
omega_x = omega(2,:);
omega_y = omega(3,:);
omega_z = omega(4,:);

num_samples = length(times);
% q_k = fromOmegaToQ([omega_x(1); omega_y(1); omega_z(1)], [dt])';
q_next_k = q_k;

for i = 1:num_samples - 1
    
    % first Runge-Kutta coefficient
    q_i_1 = q_k;
    OMEGA_omega_t_k = ...
        [0           -omega_x(i)  -omega_y(i)  -omega_z(i);
        omega_x(i)   0            omega_z(i)   -omega_y(i);
        omega_y(i)   -omega_z(i)  0            omega_x(i);
        omega_z(i)   omega_y(i)   -omega_x(i)  0          ];
    k_1 = (1/2)*OMEGA_omega_t_k*q_i_1;
    
    % second Runge-Kutta coefficient
    q_i_2 = q_k + dt(i)*(1/2)*k_1;
    OMEGA_omega_t_k_plus_half_dt = ...
        [0                                -(omega_x(i) + omega_x(i + 1))/2   -(omega_y(i) + omega_y(i + 1))/2  -(omega_z(i) + omega_z(i + 1))/2;
        (omega_x(i) + omega_x(i + 1))/2   0                                  (omega_z(i) + omega_z(i + 1))/2   -(omega_y(i) + omega_y(i + 1))/2;
        (omega_y(i) + omega_y(i + 1))/2   -(omega_z(i) + omega_z(i + 1))/2   0                                 (omega_x(i) + omega_x(i + 1))/2;
        (omega_z(i) + omega_z(i + 1))/2   (omega_y(i) + omega_y(i + 1))/2    -(omega_x(i) + omega_x(i + 1))/2  0                              ];
    k_2 = (1/2)*OMEGA_omega_t_k_plus_half_dt*q_i_2;
    
    % third Runge-Kutta coefficient
    q_i_3 = q_k + dt(i)*(1/2)*k_2;
    OMEGA_omega_t_k_plus_half_dt = ...
        [0                                -(omega_x(i) + omega_x(i + 1))/2   -(omega_y(i) + omega_y(i + 1))/2  -(omega_z(i) + omega_z(i + 1))/2;
        (omega_x(i) + omega_x(i + 1))/2   0                                  (omega_z(i) + omega_z(i + 1))/2   -(omega_y(i) + omega_y(i + 1))/2;
        (omega_y(i) + omega_y(i + 1))/2   -(omega_z(i) + omega_z(i + 1))/2   0                                 (omega_x(i) + omega_x(i + 1))/2;
        (omega_z(i) + omega_z(i + 1))/2   (omega_y(i) + omega_y(i + 1))/2    -(omega_x(i) + omega_x(i + 1))/2  0                              ];
    k_3 = (1/2)*OMEGA_omega_t_k_plus_half_dt*q_i_3;
    
    % forth Runge-Kutta coefficient
    q_i_4 = q_k + dt(i)*1*k_3;
    OMEGA_omega_t_k_plus_dt = ...
        [0               -omega_x(i + 1)  -omega_y(i + 1)  -omega_z(i + 1);
        omega_x(i + 1)   0                omega_z(i + 1)   -omega_y(i + 1);
        omega_y(i + 1)   -omega_z(i + 1)  0                omega_x(i + 1);
        omega_z(i + 1)   omega_y(i + 1)   -omega_x(i + 1)  0          ];
    k_4 = (1/2)*OMEGA_omega_t_k_plus_dt*q_i_4;
   
    q_next_k = q_k + dt(i)*((1/6)*k_1 + (1/3)*k_2 + (1/3)*k_3 + (1/6)*k_4);
    
    q_next_k = q_next_k/norm(q_next_k);

    q_k = q_next_k;
    
end
end