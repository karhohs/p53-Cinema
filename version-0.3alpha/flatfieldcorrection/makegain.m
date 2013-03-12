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
function []=makegain(chan,ffpath,dirCon_ff)
disp(['making gain image for the ' chan ' channel...'])
% identify all the exposure images
Temp=cell([1,length(dirCon_ff)]); %Initialize cell array
% ----- Identify the legitimate stacks -----
expr=[chan '(?=_\d+)'];
i=1;
for j=1:length(dirCon_ff)
    Temp2=regexp(dirCon_ff(j).name,expr,'match','once','ignorecase');
    if Temp2
        Temp{i}=dirCon_ff(j).name;
        i=i+1;
    end
end
% ----- Remove empty cells -----
Temp(i:end)=[];
%identify the length of exposure for each image
expr=[chan '_(\d+)'];
exposure = zeros(size(Temp));
flatfieldIM = cell(size(Temp));
for i=1:length(Temp)
    [~, temp_num] = regexp(Temp{i},expr,'match','once','tokens');
    exposure(i) = str2double(temp_num);
    if exposure == 0
        ind = i;
    end
    info = imfinfo([ffpath,'\',chan,'_0'],'tif');
    flatfieldIM{i}=double(imread([ffpath,'\',Temp{i}],'tif','Info',info));
end

%% weight the dark image by 5
% We have a high confidence in the dark field image and want the line to
% pass through this point more than any other. Therefore it is weighted.
Temp = zeros(1,5);
exposure = [exposure Temp];
Temp = cell(1,5);
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
        [~,b]=leastsquaresfit(x,y);
        gainIM(j,k)=b;
    end
end
gainIM=gainIM/mean(mean(gainIM));
h = fspecial('average',[15 15]);
gainIM=imfilter(gainIM,h,'replicate');
max_temp=max(max(gainIM));
max_temp=round(max_temp*1000)/1000;
im_temp=gainIM*65536/max_temp;
im_temp=uint16(im_temp);
max_temp=sprintf('%d',max_temp*1000);
imwrite(im_temp,[ffpath,'\',chan,'_gain',max_temp,'.tif'],'tif','Compression','none');
end