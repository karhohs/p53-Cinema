function [] = p53LightsCameraAction(centroid)
%centroid = a 2 column matrix with the y value in column 1 and the x value in column 2.

%Assumption #1: the source images are stored as PNGs
%Assumption #2: the output format is .mp4

%Where are the source images for the movie?
pathbf = 'G:\KWKDocuments\MATLAB\p53CinemaOutput\20120818_png\position1\wavelength1_Brightfield';
pathyfp = 'G:\KWKDocuments\MATLAB\p53CinemaOutput\20120818_png\position1\wavelength1_YFP';
pathout = 'G:\KWKDocuments\MATLAB\p53CinemaOutput';
%What is the size of the desired output?
widthOut = 1136; %iPhone5
heightOut = 640; %iPhone5

%What is the size of the input image?
widthIn = 672; %Hamamatsu binning2
heightIn = 512; %Hamamatsu binning2

%Will there be magnification of the source image?
magnification = 4; %4x is nice for binning2 input and iPhone5 output

%One piece of information needed to create a movie will be the (x,y) of the
%upperleft corner of the subimage and the width and height. It is assumed
%that the upperleft corner of the original image is the origin. Another
%piece of information are the names of the files that hold the images.

%First find the path a cell has taken and smooth it a bit.
xsmooth = round(smooth(centroid(:,2),10));
ysmooth = round(smooth(centroid(:,1),10));

%The cell will be in the center of the movie
widthout2 = widthOut/magnification;
if mod(widthout2,2) %odd
    wl=(widthout2-1)/2;
    wr=(widthout2-1)/2;
else %even
    wl=(widthout2-2)/2;
    wr=(widthout2)/2;
end
heightout2 = heightOut/magnification;
if mod(heightout2,2) %odd
    ht=(heightout2-1)/2;
    hb=(heightout2-1)/2;
else %even
    ht=(heightout2-2)/2;
    hb=(heightout2)/2;
end

xout = zeros(size(xsmooth));
yout = zeros(size(ysmooth));
xout2 = zeros(size(xsmooth));
yout2 = zeros(size(ysmooth));
for i=1:length(xsmooth)
    xout(i) = xsmooth(i)-wl;
    if xout(i)<1
        xout(i)=1;
    elseif (xsmooth(i)+wr)>widthIn
        xout(i)=widthIn-wr-wl;
    end
    xout2(i) = xout(i)+wr+wl;
    yout(i) = ysmooth(i)-ht;
    if yout(i)<1
        yout(i)=1;
    elseif (ysmooth(i)+hb)>heightIn
        yout(i)=heightIn-ht-hb;
    end
    yout2(i) = yout(i)+ht+hb;
end

%find the list of images to be used to make this movie
filenamebfproto = '20120818_s1_w1Brightfield_t1_z1.png';
filenameyfpproto = '20120818_s1_w2YFP_t1_z1.png';
lenxout = length(xout)*4-3;
indx = 1:4:lenxout;
bfObj = VideoWriter(fullfile(pathout,'s1_cell8_bf.mp4'),'MPEG-4');
bfObj.FrameRate = 15;
open(bfObj);
yfpObj = VideoWriter(fullfile(pathout,'s1_cell8_yfp.mp4'),'MPEG-4');
yfpObj.FrameRate = 15;
open(yfpObj);


for i=1:length(xout)
    filename = regexprep(filenamebfproto,'(?<=_t)\d+',num2str(indx(i)));
    I1 = imread(fullfile(pathbf,filename));
    filename = regexprep(filenameyfpproto,'(?<=_t)\d+',num2str(indx(i)));
    I2 = imread(fullfile(pathyfp,filename));
    I1 = I1(yout(i):yout2(i),xout(i):xout2(i));
    I2 = I2(yout(i):yout2(i),xout(i):xout2(i));
    [I1,I2] = rescale(I1,I2);
    writeVideo(bfObj,I1);
    writeVideo(yfpObj,I2);
    
end
close(bfObj);
close(yfpObj);

function [I1,I2] = rescale(I1,I2)
bfscale  = [2500 8000];
yfpscale = [1000, 30000];
bffactor = 255/(bfscale(2)-bfscale(1));
yfpfactor = 255/(yfpscale(2)-yfpscale(1));
for i=1:numel(I1)
    I1(i) = (I1(i)-bfscale(1))*bffactor;
    I2(i) = (I2(i)-yfpscale(1))*yfpfactor;
end
I1 = uint8(imresize(I1,4));
I2 = uint8(imresize(I2,4));

function [p53Script] = p53Screenwriter()

