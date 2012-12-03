function [] = MetamorphSlideScan2PNG(path)
% [] = my_tiffStacker(path,positions,timepoints)
% Input:
% path: a char. The path to the folder that contains the raw TIFF images from
% imageJ.
% positions (optional): an array of positive integers. The integers in this array represent the
% positions for which stacks will be created.
% timepoints (optional): a cell containing arrays of positive integers or
% an array of positive integers. If it is a cell, then it must contain an
% array for each position. If it is an array, only these timepoints will be
% taken for each position.
%
% Output:
% There is no direct argument output. Rather, stacks will be created from the
% raw images and stored in a new directory called "Stacks" at the same
% level as the path/ directory.
%
% Description:
% When an image is created in Metamorph's Multi-Dimensional-Acquisition it
% is saved in the TIFF format. An inividual image is created for each each
% stage position, wavelength, and timepoint. If z-positions are taken then
% the individual images are instead stacks of z-positions.
% When metamorph saves a .tif the filename has the following formats:
% Multi-dim. Aq. = <user defined input>_w\d+<channel name>_s\d+_t\d+.tif
% '_w' is the wavelength
% '_s' is the stage position
% '_t' is the time point (slidescan always has '_t1'; can it take time pts?)
% This function collates the images of each position and wavelength into a
% TIFF stack, because this format was necessary to view the images in
% ImageJ and enable manual segmentation of the individual cells.
%
% Other Notes:
% It is often convenient to analyze z-stacks using maximum intensity
% projection (MIP). This projects a 3D object onto a 2D plane that can then
% be analyzed with existing 2D image processing techiniques.

myConfig = openConfig();

cd([path,'\..']);
warning('off','MATLAB:tifflib:libraryWarning');
% ----- Organizing the filenames -----
%The filenames are organized within a 3D (i,j,k) matrix where each dimension
%represents either time, position, or wavelength (i=time, j=position,
%k=wavelength). This organization facilitates the collating these images
%into TIFF stacks.

FileNames=cell(2048,8); %Assumes there will not be more than 2048 positions (over a week at images every 5 minutes), or more than 8 channels. If there are, change this number. It has been pre-allocated for speed.
%FileNamesPosition,Channel)
disp(['working in ', path]); %Collect the directory contents
dirCon = dir(path);
[S,W]=deal(zeros(length(dirCon),1));
for i = 1:length(dirCon)
    if isempty(regexpi(dirCon(i).name,'_thumb')) && ...
            ~isempty(regexpi(dirCon(i).name,'(?<=(_s\d+).*)\.tif')) %Look at .tif images, but not "thumb" images
        S(i)=str2double(regexp(dirCon(i).name,'(?<=_s)\d+','match','once')); %Find the stage position number
        W(i)=str2double(regexp(dirCon(i).name,'(?<=_w)\d+','match','once')); %Find the wavelength number
        FileNames{S(i),W(i)}=dirCon(i).name; %Cleverly store the file name in a cell
    end
end
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
outpath = myConfig.output_directory.txtd8a{1};
pngpath = fullfile(outpath,sprintf('%s_png',labelText));
mkdir(pngpath);
metadatapath = fullfile(outpath,sprintf('%s_metadata',labelText));
mkdir(metadatapath);

%make the PNG images and metadata XML files

%     %Create the position directory
%     positionpath = fullfile(pngpath,sprintf('position%d',i));
%     mkdir(positionpath);
%     positionpathmeta = fullfile(metadatapath,sprintf('position%d',i));
%     mkdir(positionpathmeta);
for j = Wunique
    %Create the wavelength directory
    wavepath = fullfile(pngpath,sprintf('wavelength%d_%s',j,Wnames{j}));
    mkdir(wavepath);
    wavepathmeta = fullfile(metadatapath,sprintf('wavelength%d_%s',j,Wnames{j}));
    mkdir(wavepathmeta);
    for i = Sunique
        if ~isempty(FileNames{i,j}) %FileNames(Time,Position,Channel)
            %Load the image
            filename = fullfile(path,FileNames{i,j});
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
                
                filenamePNG = fullfile(wavepath,sprintf('%s_s%d_w%d%s_t%d_z%d.png',labelText,i,j,Wnames{j},1,h));
                imwrite(I,filenamePNG,'png','bitdepth',16);
                %Create the metadatafile
                p.filename = fullfile(wavepathmeta,sprintf('%s_s%d_w%d%s_t%d_z%d.xml',labelText,i,j,Wnames{j},1,h));
                p.labelText = labelText;
                p.zSliceText = num2str(h);
                p.positionText = num2str(i);
                p.timepointText = num2str(1);
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
                            p.stagePostitionXText = my_list.item(g).getAttribute('value').toString.toCharArray';
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


disp('cool');