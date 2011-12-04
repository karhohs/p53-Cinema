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
dirCon_ff = dir(ffpath);
%%%%%%%%%%%%%%%%%%%%%%%
% Flat Field Correction
%%%%%%%%%%%%%%%%%%%%%%%
%----- Create a new folder to hold corrected images -----
tempfoldername=regexp(stackpath,'(?<=\\)[\w ]*','match'); %Prepare to create a new folder to place background subtracted stacks
tempfoldername=[tempfoldername{end},'_ff'];
ffstackpath=[stackpath,'\..\',tempfoldername];
    mkdir(ffstackpath);
    cd(ffstackpath)
    ffstackpath=pwd;
    
%---- import directory of flatfield images ----
dirCon_ff = dir(ffpath);
%---- check for existence of correction images ----
%first, identify offset or gain images and their channel
%second, identify the blank and dark images and their channel
%identify the intersection between the sets of channels
%for the channels in the blank and dark image set, but not in the intersection, create correction images.
%---- correct stacks using correction images ----
    [flatfieldIM,fluochannels,exposure]=importflatfieldimages(ffpath); %NOTE: It is inefficient to load images and check for channel names in one function right before checking for the existence of existing offset and gain images.
    [offset,gain]=deal(cell(size(fluochannels)));
    %Use any existing flatfield images and if they don't exist create them
    for i=1:length(fluochannels)
        try
            offset{i}=double(imread([ffpath,'\',fluochannels{i},'_offset.tif']));
            flag=0; %This flag will ensure a gain file is created if no gain file exists
            for j=1:length(dirCon_ff)
                expr=['(?<=',fluochannels{i},'_gain)\d+'];
                Temp=regexp(dirCon_ff(j).name,expr,'match','once');
                if ~isempty(Temp)
                    gain{i}=double(imread([ffpath,'\',dirCon_ff(j).name]));
                    max_temp=str2double(Temp)/100;
                    gain{i}=gain{i}*max_temp/65536;
                    flag=1;
                end
            end
            if ~flag
               error('gain file does not exist');
            end
        catch ME
            disp(ME.message)
            disp(['Creating a new offset and gain image for ',fluochannels{i}])
            [offset{i},gain{i}]=findflatfield(flatfieldIM{i},exposure{i});
            im_temp=uint16(offset{i});
            imwrite(im_temp,[ffpath,'\',fluochannels{i},'_offset.tif'],'tif','Compression','none');
            max_temp=max(max(gain{i}));
            max_temp=round(max_temp*100)/100;
            im_temp=gain{i}*65536/max_temp;
            im_temp=uint16(im_temp);
            max_temp=sprintf('%d',max_temp*100);
            imwrite(im_temp,[ffpath,'\',fluochannels{i},'_gain',max_temp,'.tif'],'tif','Compression','none');
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