function []=flatfieldcorrection(stackpath,ffpath)
% [] = flatfieldcorrection(stackpath,ffpath)
% Input:
% stackpath: the location of TIFF that will be flat field corrected
% ffpath: the location of the correction images
%
% Output:
% There is no direct argument output. Rather, new stacks will be created from
% the stacks found in the stackpath directory and stored in a folder called
% [stackpath '_ff'].
%
% Description:
%
%
% Other Notes:
%
if nargin==0
    stackpath = 'C:\Users\Kyle\Documents\p53CinemaSTACKS';
    ffpath='G:\KWKDocuments\My Dropbox\p53Cinema\flatfield_20110909';
end
%----- Import stacknames -----
stacknames=importStackNames(stackpath);
%identify the channels
expr='(?<=_w\d+).*(?=_s\d+)';
Temp=cell([1,length(stacknames)]); %Initialize cell array
i=1;
for j=1:length(dirCon_stack)
    Temp2=regexp(dirCon_stack(j).name,expr,'match','once','ignorecase');
    if Temp2
        Temp{i}=Temp2;
        i=i+1;
    end
end
Temp(i:end)=[];
channels_stacks = unique(Temp); %The different channels are saved

%----- Create a new folder to hold corrected images -----
tempfoldername=regexp(stackpath,'(?<=\\)[\w ]*','match'); %Prepare to create a new folder to place background subtracted stacks
tempfoldername=[tempfoldername{end},'_ff'];
ffstackpath=[stackpath,'\..\',tempfoldername];
mkdir(ffstackpath);
cd(ffstackpath)
ffstackpath=pwd;

%Create channelTruthTable variable. The channelTruthTable variable is two
%columns with a row for each channel. The first column is the offset
%column, 1 if offset image exists 0 otherwise. The second column is the
%gain column, 1 if gain image exists 0 otherwise.
channelTruthTable = zeros(length(channels_stacks),2);
%---- check for existence of correction images ----
%import directory of flatfield images
disp(['working in ', ffpath]); %Sanity Check
dirCon_ff = dir(ffpath);
%first, identify offset  and gain images and their channel
expr='.+(?=_offset.tif)';
for j=1:length(dirCon_ff)
    Temp2=regexp(dirCon_ff(j).name,expr,'match','once','ignorecase');
    if Temp2
        for k=1:length(channels_stacks)
            if strcmpi(Temp2,channels_stacks{k})
                channelTruthTable(k,1) = 1;
            end
        end
    end
end
expr='.+(?=_gain\d+.tif)';
for j=1:length(dirCon_ff)
    Temp2=regexp(dirCon_ff(j).name,expr,'match','once','ignorecase');
    if Temp2
        for k=1:length(channels_stacks)
            if strcmpi(Temp2,channels_stacks{k})
                channelTruthTable(k,2) = 1;
            end
        end
    end
end
%create offset and gain images according to the truth table
for i=1:length(channels_stacks)
    outcome = channelTruthTable(k,1)*2 + channelTruthTable(k,2);
    switch outcome
        case 0 %No offset or gain image
            makeoffset(channels_stacks{i},ffpath);
            makegain(channels_stacks{i},ffpath);
        case 1 %just a gain image
            makeoffset(channels_stacks{i},ffpath);
        case 2 %just an offset image
            makegain(channels_stacks{i},ffpath);
        case 3 %both gain and offset exist
    end
end

%---- correct stacks using correction images ----
for j=1:length(channels_stacks)
    info = imfinfo([ffpath,'\',chan,'_offset'],'tif');
    offset = double(imread([ffpath,'\',chan,'_offset']),'tif','Info',info);
    for k=1:length(dirCon_ff)
        temp = regexp(dirCon_ff(j).name,[chan '_gain\d+'],'match','once','ignorecase');
        if ~isempty(temp)
            gainname = temp;
            break
        end
    end
    info = imfinfo([ffpath,'\',gainname],'tif');
    gain = double(imread([ffpath,'\',gainname]),'tif','Info',info);
    for i=1:length(stacknames)
        temp = regexp(stacknames{i},channels_stacks{j},'match','once');
        if ~isempty(temp)
            disp(['Flat field correction in ',stacknames{i}])
            S=readStack([stackpath,'\',stacknames{i}]);
            for k=1:length(S)
                S{k}=S{k}-offset{j};
                S{k}(S{k}<0)=0;
                S{k}=S{k}./gain{j};
            end
            Name = regexprep(stacknames{i},'(?<=_t)(\w*)(?=\.)','$1_ff');
            newStack(S,Name,ffstackpath);
        end
    end
end
end

%----- SUBFUNCTION IMPORTSTACKNAMES -----
function [Temp] = importStackNames(stackpath)
dirCon_stack = dir(stackpath);
disp(['working in ', stackpath]); %Sanity Check
expr='.*_w\d+.*_s\d+.*_t.*\.tif';
Temp=cell([1,length(dirCon_stack)]); %Initialize cell array
% ----- Identify the legitimate stacks -----
i=1;
for j=1:length(dirCon_stack)
    Temp2=regexp(dirCon_stack(j).name,expr,'match','once','ignorecase');
    if Temp2
        Temp{i}=Temp2;
        i=i+1;
    end
end
% ----- Remove empty cells -----
Temp(i:end)=[];
% for j=length(Temp):-1:1
%     if isempty(Temp{j})
%         Temp(j)=[];
%     end
% end
end

%----- SUBFUNCTION LEASTSQUARESFIT -----
function [a,b]=leastsquaresfit(x,y,j,k)
xm=mean(x);
ym=mean(y);
SSxx=sum(x.*x)-length(x)*xm^2;
SSyy=sum(y.*y)-length(y)*ym^2;
SSxy=sum(x.*y)-length(x)*xm*ym;
b=SSxy/SSxx;
a=ym-b*xm;
r2=(SSxy^2)/(SSxx*SSyy);
if r2<0.8
    disp(['not a good fit. x=', num2str(k), ' y=', num2str(j)])
end
end

%----- SUBFUNCTION XYSMOOTHEN -----
function [I]=xysmoothen(I,windowsize)
ho=size(I,1); %original height of the input image (usually 512)
wo=size(I,2); %original width of the input image (usually 672)
I2=zeros(size(I));
%ensure the window size is odd
odd=mod(windowsize,2);
if ~odd
    windowsize=windowsize+1;
end
a=1/windowsize^2;
padsize=(windowsize-1)/2;
I=padarray(I,[padsize padsize],'symmetric');
jarray=(1+padsize:1:ho+padsize);
karray=(1+padsize:1:wo+padsize);
winarray=(1:windowsize)-(padsize+1);
for j=jarray
    for k=karray
        b=0;
        for q=winarray
            for r=winarray
                b=b+I(j+q,k+r);
            end
        end
        I2(j-padsize,k-padsize)=a*b;
    end
end
I=I2;
end

%----- SUBFUNCTION MAKEOFFSET -----
function []=makeoffset(chan,ffpath)
info = imfinfo([ffpath,'\',chan,'_0'],'tif');
IM=double(imread([ffpath,'\',chan,'_0']),'tif','Info',info);
IM=xysmoothen(IM,9);
IM=floor(IM);
IM=uint16(IM);
imwrite(IM,[ffpath,'\',chan,'_offset'],'tif','Compression','none');
end

%----- SUBFUNCTION MAKEGAIN -----
function []=makegain(chan,ffpath,dirCon_ff)
disp(['making gain image for the ' chan 'channel...'])
%identify all the exposure images
Temp=cell([1,length(dirCon_ff)]); %Initialize cell array
% ----- Identify the legitimate stacks -----
expr=[chan '(?=_\d+)'];
i=1;
for j=1:length(dirCon_ff)
    Temp2=regexp(dirCon_ff(j).name,expr,'match','once','ignorecase');
    if Temp2
        Temp{i}=dirCon_ff(j);
        i=i+1;
    end
end
% ----- Remove empty cells -----
Temp(i:end)=[];
%identify the length of exposure for each image
exposure = zeros(size(Temp));
flatfieldIM = cell(size(Temp));
for i=1:length(Temp)
    exposure(i) = str2double(regexp(Temp{i},'\d+','match','once'));
    if exposure == 0
        ind = i;
    end
    info = imfinfo([ffpath,'\',chan,'_0'],'tif');
    flatfieldIM{i}=double(imread([ffpath,'\',Temp{i}]),'tif','Info',info);
end

%weight the dark image by 5
flatfieldIM{end+1:end+5} = flatfieldIM{ind};
exposure(end+1:end+5) = 0;

%calculate the gain image
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
gainIM=xysmoothen(gainIM,5);
max_temp=max(max(gainIM));
max_temp=round(max_temp*1000)/1000;
im_temp=gainIM*65536/max_temp;
im_temp=uint16(im_temp);
max_temp=sprintf('%d',max_temp*1000);
imwrite(im_temp,[ffpath,'\',chan,'_gain',max_temp],'tif','Compression','none');
end