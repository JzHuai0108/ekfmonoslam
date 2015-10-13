function testvlfeatsift()
 dir='\\File.ceegs.ohio-state.edu\SPIN\UAV\August 08, 2013\Nikon D800\109ND800\';
fn1=[dir,'DSC_1023.jpg'];
fn2=[dir,'DSC_1024.jpg'];
im1=imread(fn1);
im2=imread(fn2);
Nimg1 = gaussdownsample(im1,4,[5 5],1.6);
Nimg2 = gaussdownsample(im2,4,[5 5],1.6);

figure; imshow(Nimg1);

sift_mosaic(Nimg1, Nimg2);  

function Nimg = gaussdownsample(img,levels,gsize,sigma);
%
%  Function to generate a gaussian-pyramid for the given input image
% image is resize to width/2^(levels-1)*height/2^(levels-1)
% Input:  
%        img: input image-matrix grayscale
%        levels:  number of levels of the pyramid, we only return the level
%        of the least size
%        gsize: size of the gaussian kernel [w h] ([5 5] normally provides a smooth output)
%        sigma:  sigma for gaussian kernel 
% Output:
%        Nimg:  the smoothed image of the smallest size
% Usage:
%      im = imread('cameraman.tif');
%      Nimg = Gscale(im,3,[5 5],1.6);
%      figure; imshow(Nimg);
Nimg=img;
if(levels==1)
    return;
end
%guassian filter  with a sigma=1.6
g = fspecial('gaussian',gsize,sigma);
%pyramid
for i = 2:levels
        %perform guassian filtering
        im = imfilter(Nimg,g,'conv');
        %perform downsampling (horizontal)
        im1 = im(:,1:2:size(Nimg,2),:);
        %vertical
        im2 = im1(1:2:size(Nimg,1),:,:);
        %store it in a struct format
        Nimg = im2;    
end
        
%End
