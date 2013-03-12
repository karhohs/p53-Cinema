%% flatfieldcorrection
% Using dark field and flat field images collected prior to or after the
% time-lapse fluorescent images, biases created by a spatial heterogenity
% in light source intensity and light path absorbance, and the bias added
% to the image by the camera during acquisition, are removed in a precise,
% quantitative fashion.
%
%   [] = flatfieldcorrection(path,ffpath)
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
% The light measured by the camera is thought to be proportionally related
% to the length of time exposed to the incident light. Then there is an
% offset created by the measurement device that must also be accounted for.
% This describes a linear relationship.
%
% $$y = \beta t + C$$
%
% Here, $y$ is the pixel measurement; $\beta$ is a combination of flux,
% photons per pixel area per time, and light path absorbance, unitless
% measure between 0 an 1; $t$ is time; $C$ is the dark field offset. The
% two parameters $\beta$ and $C$ can be found by fitting a line for each
% pixel using images taken at several different exposures.
%
% Other Notes:
% Assumes the image files have been converted to the p53Cinema PNG format
function []=flatfieldcorrection(path,ffpath,varargin)

p = inputParser;
p.addRequired('path', @(x)ischar(x));
p.addRequired('ffpath', @(x)ischar(x));
p.parse(path,ffpath, varargin{:});

%----- Import PNG paths -----
%Check whether the computer is a mac
if ismac
    error('ffcor:notApc','Flatfield correction is currently configured to run on a pc');
    %[~,filepaths]=system('ls -Rp */*.png');
elseif ispc
    %the 'dir' command in Matlab does not search subfolders. However,
    %there is a MS-DOS command that will do the trick!
    %'/B' uses bare format (no heading information or summary)
    %'/S' displays files in specified directory and all subdirectories
    cd(path)
    [~,filepaths] = system('dir /S /B *.png');
    filepaths = textscan(filepaths,'%s');
    filepaths = filepaths{1};
else
    error('ffcor:notApc','Flatfield correction is currently configured to run on a pc');
end
disp(['working in ', path]); %Sanity Check

%import directory of flatfield images
disp(['working in ', ffpath]); %Sanity Check
dirCon_ff = dir(ffpath);
%only the channels with flat field images can be corrected. identify those
%channels
expr='.+(?=_)';
Temp=cell([1,length(dirCon_ff)]); %Initialize cell array
i=1;
for j=1:length(dirCon_ff)
    Temp2=regexp(dirCon_ff(j).name,expr,'match','once','ignorecase');
    if Temp2
        Temp{i}=Temp2;
        i=i+1;
    end
end
Temp(i:end)=[];
channels_stacks = unique(Temp); %The different channels are saved

if isempty(channels_stacks)
    error('fltfldcrct:noFF','There are no matching flatfield images.')
end
%----- Create a new folder to hold corrected images -----
foldernameIN=regexp(path,'(?<=\\)[\w ]*','match'); %Prepare to create a new folder to place background subtracted stacks
foldernameIN = foldernameIN{end};
foldernameOUT=[foldernameIN,'_ff'];
ffstackpath=regexprep(stackpath,foldernameIN,foldernameOUT);
mkdir(ffstackpath);

%Create channelTruthTable variable. The channelTruthTable variable is two
%columns with a row for each channel. The first column is the offset
%column, 1 if offset image exists 0 otherwise. The second column is the
%gain column, 1 if gain image exists 0 otherwise.
channelTruthTable = zeros(length(channels_stacks),2);
%---- check for existence of correction images ----
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
    outcome = channelTruthTable(i,1)*2 + channelTruthTable(i,2);
    switch outcome
        case 0 %No offset or gain image
            makeoffset(channels_stacks{i},ffpath);
            makegain(channels_stacks{i},ffpath,dirCon_ff);
        case 1 %just a gain image
            makeoffset(channels_stacks{i},ffpath);
        case 2 %just an offset image
            makegain(channels_stacks{i},ffpath,dirCon_ff);
        case 3 %both gain and offset exist
    end
end

%---- correct stacks using correction images ----
dirCon_ff = dir(ffpath);
for j=1:length(channels_stacks)
    info = imfinfo([ffpath,'\',channels_stacks{j},'_offset'],'tif');
    offset = double(imread([ffpath,'\',channels_stacks{j},'_offset'],'tif','Info',info));
    offset = scale12to16bit(offset);
    for k=1:length(dirCon_ff)
        temp = regexp(dirCon_ff(k).name,[channels_stacks{j} '_gain\d+'],'match','once','ignorecase');
        if ~isempty(temp)
            gainname = temp;
            break
        end
    end
    info = imfinfo([ffpath,'\',gainname],'tif');
    gain = double(imread([ffpath,'\',gainname],'tif','Info',info));
    expr='(?<=_gain)\d+';
    max_temp=regexp(gainname,expr,'match','once');
    max_temp=str2double(max_temp)/1000;
    gain=gain*max_temp/65536;
    
    for i=1:length(filepaths)
        [~,pngname,~] = fileparts(filepaths{i});
        temp = regexp(pngname,channels_stacks{j},'match','once','ignorecase');
        if ~isempty(temp)
            disp(['Flatfield correcting ',filepaths{i}])
            IM = double(imread(filepaths{i}));
            Name = regexprep(filepaths{i},foldernameIN,foldernameOUT);
            [pngpath,~,~] = fileparts(Name);
            %verify that the input image is the correct size, otherwise no
            %output image will be created.
            if size(IM)~=size(offset)
                continue
            end
            IM = IM-offset;
            IM(IM<0) = 0;
            IM = IM./gain;
            IM = uint16(IM);
            %check for the existence of the folder where the new image will be placed.
            if ~isdir(pngpath)
                mkdir(pngpath);
            end
            %Create the PNG
            imwrite(IM,Name,'png','bitdepth',16);
        end
        %update metadata?
        
    end
end
%% Create or append a log file
% 
fid = fopen(fullfile(path,'log.txt'),'a+');
fprintf(fid,'%s: flatfieldcorrection: PNG images from the folder ''%s'' were flatfield corrected.\r\n\r\n',date,path);
end




