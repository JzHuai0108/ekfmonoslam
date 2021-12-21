addpath('..\mexopencv\');
dir='G:\data\20130808\109ND800_080813\';
prefix=sprintf('%sDSC_', dir);

img0 = takeImage_v002(prefix, 1023, 'jpg', 3);
mask=uint8(ones(size(img0,1),size(img0,2)));
mask(:,:)=1;
excluded_band=100;
mask(1:excluded_band,:)=0;
mask(end-excluded_band+1:end,:)=0;
mask(:,1:excluded_band)=0;
mask(:,end-excluded_band+1:end)=0;

 prevPts = cv.goodFeaturesToTrack(img0,'MaxCorners', 200, ...
 'QualityLevel', 0.02,'MinDistance', 10,'Mask', mask, 'BlockSize', 3, 'UseHarrisDetector', 0, 'K',0.04);
prevPts=cv.cornerSubPix(img0, prevPts);

% the argument order and vacancy does not matter in matlab in constrast to
% opencv C++
corners=zeros(2,size(prevPts,2));
    for jerk=1:size(prevPts,2)
        if(~isempty(prevPts{jerk}))
   corners(:,jerk)= prevPts{jerk}';  
        end
    end

% increase lambda apparently increase number of corners
%larger gaussian window sz2 will slightly decrease the number of corners

figure(1); clf;
imshow(img0);
hold on;            %# Add subsequent plots to the image
%plot(corners(2,:),corners(1,:),'go');  %# NOTE: x_p and y_p are switched (see note below)!
plot(corners(1,:),corners(2,:),'ro','MarkerSize',5);

hold off;
[m,n,p]=size(img0);
winsz1=7;
sigma1=10;

imglast=img0(:,:,1);
for i=1024:1055
    im=takeImage_v002(prefix,i, 'jpg', 3);
            maxedge=max(size(im));
                while(maxedge>1000)
                    im=impyramid(im, 'reduce');% 8UC1
                    maxedge=round(maxedge/2);
                end
%     img= double((read(myVid,i-1)));
%     img2=double((read(myVid,i)));
%     for band=1:size(img,3)
%         img(:,:,band)=simgauss(img(:,:,band), winsz1, sigma1);
%         img2(:,:,band)=simgauss(img2(:,:,band), winsz1, sigma1);
%     end   
    [prevPts, status, err] = cv.calcOpticalFlowPyrLK(imglast, im, prevPts, 'WinSize', [23, 23], ...
        'MaxLevel', 2,'Criteria', struct('type', 'Count+EPS', 'maxCount', 20, 'epsilon', 0.01));
    for jerk=1:length(status)
        if(status~=0)
        end
    end
    corners=zeros(2,size(prevPts,2));
    for jerk=1:size(prevPts,2)
        if(~isempty(prevPts{jerk}))
   corners(:,jerk)= prevPts{jerk}';  
        end
    end
    
   % corners=SimpKLT(img,img2,corners,3); 
    imshow(im);
    hold on
    plot(corners(1,:),corners(2,:),'ro','MarkerSize',5);
    hold off
    drawnow;
   imglast=im;
  
end

% close(writerObj);
