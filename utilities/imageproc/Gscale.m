function Nimg = Gscale(img,levels,gsize,sigma);
%
%  Function to generate a gaussian-pyramid for the given input image
%
% Input:  
%        img: input image-matrix grayscale
%        levels:  number of levels of the pyramid 
%        gsize: size of the gaussian kernel [w h] ([5 5] normally provides a smooth output)
%        sigma:  sigma for gaussian kernel 
% Output:
%        Nimg:  is a struct consisting of images from each level
%            :  Nimg.img;
% Usage:
%      im = imread('cameraman.tif');
%      Nimg = Gscale(im,3,[5 5],1.6);
%      i = 2; %select a level
%      figure; imshow(Nimg(i).img);
%
% Author: Pranam Janney                     Date: 24th July 2006 15:39
% Email: pranamjanney@yahoo.com
%
% Revised Version 1.0.1                     Date: 04th August 2006, 10:50
%


%guassian filter  with a sigma=1.6
g = fspecial('gaussian',gsize,sigma);

%pyramid
for i = 1:levels
    if i == 1
        im = imfilter(img,g,'conv');
        Nimg(i).img = im;
    else 
        %perform guassian filtering
        im = imfilter(Nimg(i-1).img,g,'conv');
        %perform downsampling (horizontal)
        im1 = im(:,1:2:size(Nimg(i-1).img,2));
        %vertical
        im2 = im1(1:2:size(Nimg(i-1).img,1),:);
        %store it in a struct format
        Nimg(i).img = im2;
    end
end
        
%End
        
