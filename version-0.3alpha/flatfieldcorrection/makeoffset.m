%% makeoffset
% The offset image is an approximation of the dark field. The dark field
% image is the signal captured by the camera when there is no light on the
% sensor. This function will smooth the image captured with a zero
% exposure.
%
%   [] = makeoffset(path,ffpath)
%
%%% Input:
% * chan: the name of the channel to have an offset image created.
% * ffpath: the location of the correction images
%
%%% Output:
% An offset image is created and written to the directory that contains the
% flat field images.
%
%%% Description:
%
%
% Other Notes:
% 
function []=makeoffset(chan,ffpath)
disp(['making offset image for the ' chan ' channel...'])
fname = fullfile(ffpath,strcat(chan,'_0'));
info = imfinfo(fname,'tif');
IM=double(imread(fname,'tif','Info',info));
%%
% scale from 12-bit image to 16-bit image
IM = uint16(IM);
IM = bitshift(IM,4);
%% smooth the image
% Images from the lab typically end up being 1344 x 1024 or 672 x 512,
% depending on whether or not there is binning. The size of the image will
% influence the size of the filters used to smooth the image.
if info.Width == 1344
    IM = medfilt2(IM,[31,31],'symmetric'); %median filters are good for salt and pepper noise like that seen in the darkfield image
    h = fspecial('average',[31 31]); %the average filter removes abrupt spatial changes in intensity
    IM=imfilter(IM,h,'replicate'); 
else
    IM = medfilt2(IM,[15,15],'symmetric');
    h = fspecial('average',[15 15]);
    IM=imfilter(IM,h,'replicate');
end
imwrite(IM,fullfile(ffpath,strcat(chan,'_offset.tif')),'tif','Compression','none');
end