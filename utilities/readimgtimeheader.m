function [fimgtime, imgepoch, initIm]=readimgtimeheader(imgtimefile,preimutime, initIm)
fimgtime=fopen(imgtimefile,'rt');
if(fimgtime~=-1)
    % Discard header  
    fgetl(fimgtime);  
    hstream= fgetl(fimgtime);  
    camel=sscanf(hstream,'%f,%f'); 
    % find the data that has greater index and greater time epoch
    while (camel(2)<=preimutime||camel(1)<initIm)
        hstream= fgetl(fimgtime);
        if(hstream==-1)
            error('None image exists with large image epoch and frame Id!');
        end
        camel=sscanf(hstream,'%f,%f'); 
    end
    imgepoch=camel(2);
    initIm=camel(1);
else imgepoch=inf;
end
