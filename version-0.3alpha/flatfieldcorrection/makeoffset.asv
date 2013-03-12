%% makeoffset
% The offset image is an approximation of the dark field. The dark field
% image is the signal captured by the camera when there is no light on the
% sensor. This function will smooth the image captured with a zero
% exposure.
%
%   [] = makeoffset(path,ffpath)
%
%%% Input:
% * path: the location of outpath created by |importFromMetamorphMDA|. The
% PNG images within will be flat field corrected. The function
% |importFromMetamorphMDA| should be run beforehand.
% * ffpath: the location of the correction images
%
%%% Output:
% There is no direct argument output. Existing PNG images in the _path_ are
% flatfield corrected.
%
%%% Description:
%
%
% Other Notes:
% Assumes the image files have been converted to the PNG format
function []=makeoffset(chan,ffpath)
disp(['making offset image for the ' chan ' channel...'])
info = imfinfo([ffpath,'\',chan,'_0'],'tif');
IM=double(imread([ffpath,'\',chan,'_0'],'tif','Info',info));
h = fspecial('gaussian',[21 21],3.5);
IM=imfilter(IM,h,'replicate');
IM=uint16(IM);
imwrite(IM,[ffpath,'\',chan,'_offset.tif'],'tif','Compression','none');
end