%% MetamorphMDA2PNG
% Converts a folder of TIFF images created via Metamorph into PNG images.
% An XML file containing metadata about all the images is also created.
%
% PNG images are preferred to TIFF images, because the format is simpler
% and more standardized. Also, the PNG has better lossless compression.
%
%   [] = my_tiffStacker(path,positions,timepoints)
%
%%% Input
% * path: a char. The path to the folder that contains the orignal TIFF
% images from metamorph.
%
%%% Output:
% There is no direct argument output. Rather, a new folder will be created
% that contains PNG images and other files that will later be accessed by
% p53Cinema.
%
%%% Detailed Description
% When an image is created in Metamorph's Multi-Dimensional-Acquisition it
% is saved in the TIFF format. An inividual image is created for each each
% stage position, wavelength, and timepoint. If z-positions are taken then
% the individual images are instead stacks of z-positions.
% When metamorph saves a .tif the filename has the following formats:
%
% |"user defined input"_w\d+"channel name"_s\d+_t\d+.tif|
%
% * '_w' is the wavelength
% * '_s' is the stage position
% * '_t' is the time point
%
%%% Other Notes
% 
function [] = MetamorphMDA2PNG(path)
cd([path,'\..']);
outpath = fullfile(pwd,'p53Cinema');
warning('off','MATLAB:tifflib:libraryWarning');
%% Organizing the filenames
% The filenames are organized within a 3D (i,j,k) matrix where each
% dimension represents either time, position, or wavelength (i=time,
% j=position, k=wavelength). This organization facilitates the collating
% these images into TIFF stacks.
%
%   FileNames(Time,Position,Channel);
FileNames=cell(2048,256,8); %Assumes there will not be more than 2048 timepoints (over a week at images every 5 minutes), 256 positions, or more than 8 channels. If there are, change this number. It has been pre-allocated for speed.
%% 
% Collect the directory contents
dirCon = dir(path);
[S,T,W]=deal(zeros(length(dirCon),1));
for i = 1:length(dirCon)
    if isempty(regexpi(dirCon(i).name,'_thumb')) && ...
            ~isempty(regexpi(dirCon(i).name,'(?<=(_s\d+).*)\.tif'))%Look at .tif images, but not "thumb" images
        T(i)=str2double(regexp(dirCon(i).name,'(?<=_t)\d+','match','once')); %Find the time point number
        S(i)=str2double(regexp(dirCon(i).name,'(?<=_s)\d+','match','once')); %Find the stage position number
        W(i)=str2double(regexp(dirCon(i).name,'(?<=_w)\d+','match','once')); %Find the wavelength number
        FileNames{T(i),S(i),W(i)}=dirCon(i).name; %Cleverly store the file name in a cell
    end
end
Tunique = unique(T)';
Tunique = Tunique(Tunique~=0);
Sunique = unique(S)';
Sunique = Sunique(Sunique~=0);
Wunique = unique(W)';
Wunique = Wunique(Wunique~=0);

%find the names of the wavelengths
Wnames = cell(size(Wunique));
for i=1:length(Wunique)
    ind_findwavelength = find(W==Wunique(i),1,'first');
    filename_findwavelength = dirCon(ind_findwavelength).name;
    t = Tiff(fullfile(path,filename_findwavelength),'r');
    metadata = t.getTag('ImageDescription');
    t.close;
    fid = fopen('t3mp.xml','w');
    fprintf(fid,'%s',metadata);
    fclose(fid);
    xdoc = xmlread('t3mp.xml');
    node_root = xdoc.getDocumentElement;
    my_list = node_root.getElementsByTagName('prop');
    for j=0:my_list.getLength-1
        if my_list.item(j).hasAttributes
            if strcmp(my_list.item(j).getAttribute('id').toString.toCharArray','image-name')
                Wnames{i} = my_list.item(j).getAttribute('value').toString.toCharArray';
                %edit the name to remove unwanted characters and words
                Wnames{i} = regexprep(Wnames{i},'\s',''); %Remove all not(alphabetic, numeric, or underscore) characters
                Wnames{i} = regexprep(Wnames{i},'tocamera','','ignorecase'); %remove 'tocamera' if present b/c it is not informative
                Wnames{i} = regexprep(Wnames{i},'camera','','ignorecase'); %remove 'camera' if present b/c it is not informative
            end
        end
    end
    delete('t3mp.xml');
end
%find the label
labelText = regexp(dirCon(ind_findwavelength).name,'.+(?=_w)','match','once');
%create the root directory that will store the PNG images and metadata
pngpath = fullfile(outpath,'png');
mkdir(pngpath);
%% make the PNG images
%
for i = Sunique
    %Create the position directory
    positionpath = fullfile(pngpath,sprintf('position%d',i));
    mkdir(positionpath);
    positionpathmeta = fullfile(metadatapath,sprintf('position%d',i));
    mkdir(positionpathmeta);
    for j = Wunique
        %Create the wavelength directory
        wavepath = fullfile(positionpath,sprintf('wavelength%d_%s',j,Wnames{j}));
        mkdir(wavepath);
        wavepathmeta = fullfile(positionpathmeta,sprintf('wavelength%d_%s',j,Wnames{j}));
        mkdir(wavepathmeta);
        for k = Tunique
            if ~isempty(FileNames{k,i,j}) %FileNames(Time,Position,Channel)
                %Load the image
                filename = fullfile(path,FileNames{k,i,j});
                t = Tiff(filename,'r');
                %Count how many images are stored in the TIFF file. This is
                %important because there may be >1 z-slice.
                numberOfTIFFdir = 1;
                numberOfTIFFflag = true;
                while numberOfTIFFflag
                    if t.lastDirectory
                        t.setDirectory(1);
                        numberOfTIFFflag = false;
                    else
                        t.nextDirectory;
                        numberOfTIFFdir = numberOfTIFFdir + 1;
                    end
                end
                for h = 1:numberOfTIFFdir
                    %Create the PNG
                    %Load the image
                    I = t.read;
                    
                    %The Hamamatsu cameras for the closet scope and curtain scope create images
                    %with 12-bit dynamic range. However, the TIFF format that stores these
                    %images uses a 16-bit format. Viewing a 12-bit image in a 16-bit format on
                    %a computer monitor is compromised by the scaling being done at 16-bit. To
                    %make viewing images from the microscope easier on a computer monitor,
                    %without any compression or loss of data, the 12-bit data is shifted left
                    %4-bits to become 16-bit data. In addition, more information is kept following image processing.
                    numType = class(I);
                    switch numType
                        case 'double'
                            I = uint16(I);
                            I = bitshift(I,4);
                            I = double(I);
                        case 'uint16'
                            I = bitshift(I,4);
                    end
                    
                    filenamePNG = fullfile(wavepath,sprintf('%s_s%d_w%d%s_t%d_z%d.png',labelText,i,j,Wnames{j},k,h));
                    imwrite(I,filenamePNG,'png','bitdepth',16);
                    %Create the metadatafile
                    p.filename = fullfile(wavepathmeta,sprintf('%s_s%d_w%d%s_t%d_z%d.xml',labelText,i,j,Wnames{j},k,h));
                    p.labelText = labelText;
                    p.zSliceText = num2str(h);
                    p.positionText = num2str(i);
                    p.timepointText = num2str(k);
                    p.wavelengthTypeText = Wnames{j};
                    %Explore the image description in the TIFF
                    metadata = t.getTag('ImageDescription');
                    fid = fopen('t3mp.xml','w');
                    fprintf(fid,'%s',metadata);
                    fclose(fid);
                    xdoc = xmlread('t3mp.xml');
                    node_root = xdoc.getDocumentElement;
                    my_list = node_root.getElementsByTagName('prop');
                    for g=0:my_list.getLength-1
                        if my_list.item(g).hasAttributes
                            if strcmp(my_list.item(g).getAttribute('id').toString.toCharArray','acquisition-time-local')
                                p.timeOfAcquisitionText = my_list.item(g).getAttribute('value').toString.toCharArray';
                            elseif strcmp(my_list.item(g).getAttribute('id').toString.toCharArray','stage-label')
                                p.stageLabelText = my_list.item(g).getAttribute('value').toString.toCharArray';
                            elseif strcmp(my_list.item(g).getAttribute('id').toString.toCharArray','stage-position-x')
                                p.stagePositionXText = my_list.item(g).getAttribute('value').toString.toCharArray';
                            elseif strcmp(my_list.item(g).getAttribute('id').toString.toCharArray','stage-position-y')
                                p.stagePositionYText = my_list.item(g).getAttribute('value').toString.toCharArray';
                            end
                        end
                    end
                    delete('t3mp.xml');
                    createMetadata(p);
                    if ~t.lastDirectory
                        t.nextDirectory;
                    end
                end
                t.close;
            end
        end
    end
end