%% makegain
% The gain image contains weights to counter the natural bias found in
% every pixel of the image. It is assumed that the measurement of light at
% each pixel is related to the true incident light proportionally through
% an absorbance constant.
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
function []=makegain(chan,ffpath)
disp(['making gain image for the ' chan ' channel...'])
dirCon_ff = dir(ffpath);
% identify all the exposure images
Temp=cell([1,length(dirCon_ff)]); %Initialize cell array
% ----- Identify the legitimate stacks -----
expr=strcat(chan,'(?=_\d+)');
i=1;
for j=1:length(dirCon_ff)
    Temp2=regexp(dirCon_ff(j).name,expr,'match','once','ignorecase');
    Temp3 = regexp(dirCon_ff(j).name,strcat(chan,'(?=_offset)'),'match','once','ignorecase');
    if Temp2
        Temp{i}=dirCon_ff(j).name;
        i=i+1;
    end
    if Temp3
        Temp{i}=dirCon_ff(j).name;
        ind = i;
        i=i+1;
    end
end
% ----- Remove empty cells -----
Temp(i:end)=[];
%identify the length of exposure for each image
expr=strcat(chan,'_(\d+)');
exposure = zeros(size(Temp));
flatfieldIM = cell(size(Temp));
for i=1:length(Temp)
    if i~=ind
        [~, temp_num] = regexp(Temp{i},expr,'match','once','tokens');
        exposure(i) = str2double(temp_num);
    end
    info = imfinfo(fullfile(ffpath,Temp{i}),'tif');
    flatfieldIM{i}=double(imread(fullfile(ffpath,Temp{i}),'tif','Info',info));
end
flatfieldIM{ind} = flatfieldIM{ind}/16; %to compensate for the offset image being shifted over by 4 bits.
%% weight the dark image
% We have a high confidence in the dark field image and want the line to
% pass through this point more than any other. Therefore it is weighted.
Temp = zeros(1,4);
exposure = [exposure Temp];
Temp = cell(1,4);
flatfieldIM = [flatfieldIM Temp];
for i=0:4
    flatfieldIM{end-i} = flatfieldIM{ind};
end
%% calculate the gain image
[hei,wid]=size(flatfieldIM{1});
gainIM=zeros(size(flatfieldIM{1}));
for j=1:hei
    for k=1:wid
        [x,y]=deal(zeros(length(exposure),1));
        for i=1:length(exposure)
            y(i)=flatfieldIM{i}(j,k);
            x(i)=exposure(i);
        end
        [~,b]=leastsquaresfit(x,y,j,k);
        gainIM(j,k)=b;
    end
end
gainIM=gainIM/mean(mean(gainIM));
%% smooth the image
% Images from the lab typically end up being 1344 x 1024 or 672 x 512,
% depending on whether or not there is binning. The size of the image will
% influence the size of the filters used to smooth the image.
if info.Width == 1344
    h = fspecial('average',[31 31]);
    gainIM=imfilter(gainIM,h,'replicate');
else
    h = fspecial('average',[15 15]);
    gainIM=imfilter(gainIM,h,'replicate');
end
%%
% The image is normalized by the mean. But numbers between 0 and 1 cannot
% be directly stored in an 16-bit image. Therefore, the weights are scaled
% and that scaling factor is saved in the filename, so that it can be
% inverted later on.
max_temp=max(max(gainIM));
max_temp=round(max_temp*1000)/1000;
im_temp=gainIM*65536/max_temp;
im_temp=uint16(im_temp);
max_temp=sprintf('%d',max_temp*1000);
imwrite(im_temp,fullfile(ffpath,strcat(chan,'_gain',max_temp,'.tif')),'tif','Compression','none');
end