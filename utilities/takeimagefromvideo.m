% input: the handle to the video, the frmIndex 1 based, the downscale for
% the image, each edge is reduced to 1/(1<<scale)
function im=takeimagefromvideo(myVid, frmindex, scale)
if(nargin<3)
    scale=0;
end
imRGB= read(myVid,frmindex);
im=imRGB(:,:,1);
loop=0;
while(loop<scale)
    im=impyramid(im, 'reduce');
    loop=loop+1;  
end