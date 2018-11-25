addpath('..\mexopencv\');
addpath('..\utilities');
sequencePath ='G:\data\20130808\casio2_3430.MOV';
myVid= VideoReader(sequencePath);
img0 = takeimagefromvideo(myVid, 12180, 1);

mask=uint8(ones(size(img0,1),size(img0,2)));
mask(:,:)=1;
excluded_band=40;
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
for i=12181:13000
    im = takeimagefromvideo(myVid, i, 1);
    [prevPts, status, err] = cv.calcOpticalFlowPyrLK(imglast, im, prevPts, 'WinSize', [23, 23], ...
        'MaxLevel', 2,'Criteria', struct('type', 'Count+EPS', 'maxCount', 20, 'epsilon', 0.01));
    corners=zeros(2,size(prevPts,2));
    zane=0;
    for jerk=1:size(prevPts,2)
        if(status(jerk)==1&&~isempty(prevPts{jerk}))
            zane=zane+1;
            corners(:,zane)= prevPts{jerk}';
        end
    end
    if(zane<size(prevPts,2)/2)% add more points if lose too much
        prevPts = cv.goodFeaturesToTrack(img0,'MaxCorners', 200, ...
    'QualityLevel', 0.02,'MinDistance', 10,'Mask', mask, 'BlockSize', 3, 'UseHarrisDetector', 0, 'K',0.04);
    prevPts=cv.cornerSubPix(img0, prevPts);
    end
    % corners=SimpKLT(img,img2,corners,3);
    imshow(im);
    hold on
    plot(corners(1,1:zane),corners(2,1:zane),'ro','MarkerSize',5);
    hold off
    drawnow;
    imglast=im;    
end
