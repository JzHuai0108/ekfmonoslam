% draw the frame id (1 based ) and GPS TOW timestamp onto the image
function drawTimeOntoVideo()
close all;
resdir='..\20130808\';
sequencePath =[resdir, 'casio2_3430.MOV'];
imgtimefile=[resdir, 'casio2timestamps.txt'];
initIm =12464; % 415530 GTOW
endIm=12986;
[fimgtime, imgepoch, initIm]=readimgtimeheader(imgtimefile,0, initIm);
hstream=sprintf('%d,%.5f', initIm, imgepoch);
myVid= VideoReader(sequencePath);
figure(1);
tic;
for x = initIm : endIm
    imRGB=read(myVid,x); %78612;
    imwrite(imRGB,strcat('F:\relaylatest\20130808\garage\casio2_',num2str(x),'.jpg'));
    %     if(isempty(imRGB)||imgepoch>416033.2)
    %         break;
    %     else
    %display frmId and timestamp on the image
    %         text(xlbl(:), ylbl(:), lbl(:),'color','w',...
    %     'HorizontalAlignment','center','VerticalAlignment','middle');
    imshow(imRGB);
    text(50, 50, hstream);
    pause(0.02);
    % read in next image epoch
    hstream= fgetl(fimgtime);
    %         if(hstream==-1)
    %             imgepoch=inf;
    %         end
end
toc;
end