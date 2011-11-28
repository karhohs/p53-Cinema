function []=my_tiffStacker(path,positions,timepoints)
% S = my_parseXML(filename)
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
% Slide Scan = <user defined input>_w\d+_s\d+_t\d+.tif
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

cd([path,'\..']);
path2=pwd;

% ----- Organizing the filenames -----
%The filenames are organized within a 3D (i,j,k) matrix where each dimension
%represents either time, position, or wavelength (i=time, j=position,
%k=wavelength). This organization facilitates the collating these images
%into TIFF stacks.

FileNames=cell(1024,256,8); %Assumes there will not be more than 1024 timepoints, 256 positions, or more than 8 channels. If there are, change this number. It has been pre-allocated for speed.
%FileNames(Time,Position,Channel)
disp(['working in ', path]); %Collect the directory contents
dirCon = dir(path);
[Smax,Tmax,Wmax]=deal(0);
for i = 1:length(dirCon)
    if isempty(regexpi(dirCon(i).name,'_thumb')) && ...
            ~isempty(regexpi(dirCon(i).name,'(?<=(_s\d+).*)\.tif'))%Look at .tif images, but not "thumb" images
        NTP=str2double(regexp(dirCon(i).name,'(?<=_t)\d+','match','once')); %Find the time point number
        NPos=str2double(regexp(dirCon(i).name,'(?<=_s)\d+','match','once')); %Find the stage position number
        NW=str2double(regexp(dirCon(i).name,'(?<=_w)\d+','match','once')); %Find the wavelength number
        FileNames{NTP,NPos,NW}=dirCon(i).name; %Cleverly store the file name in a cell
        if NPos>Smax %The maximum number of stage positions is found
            Smax=NPos;
        end
        if NTP>Tmax %The maximum number of time points is found
            Tmax=NTP;
        end
        if NW>Wmax %The maximum number of wavelengths is found
            Wmax=NW;
        end
    end
end

% ----- process optional arguments -----
if ~exist('positions','var')
    positions = 1:Smax;
elseif min(size(positions)) > 1 || max(size(positions)) > Smax
    error('the positions vector is an invalid size')
end
if ~exist('timepoints','var')
    timepoints = cell(1,Smax);
    temp = 1:Tmax;
    for i=1:Smax
        timepoints{i} = temp;
    end
    clear('temp')
elseif isnumeric(timepoints) && min(size(timepoints)) == 1
    temp = sort(timepoints);
    timepoints = cell(1,Smax);
    for i=1:Smax
        timepoints{i} = temp;
    end
elseif min(size(timepoints)) ~= 1 || max(size(timepoints)) ~= Smax
    error('the timepoints cell vector is an invalid size')
end
mkdir(path2,'Stacks') %a new directory is generated for the stacks

%for each wave length and position, the tiff files from each time point is
%loaded and appendend to a new (multi-layer) tif file named with '_tSTACK'

% ----- Create stacks -----
timelabels = cell(Tmax,Smax,Wmax); %This cell will contain the times images were acquired
newnames = cell(Smax,Wmax);
for i=1:Wmax
    for j=positions
        for k=timepoints{j}
            if ~isempty(FileNames{k,j,i})
                %Loading the images is the most time intensive part of the
                %code
                %When metamorph saves a .tiff the filename has the following formats:
                %Multi-dim. Aq. = <user defined input>_w\d+<channel name>_s\d+_t\d+.tif
                %Slide Scan = <user defined input>_w\d+_s\d+_t\d+.tif
                %'_w' is the wavelength
                %'_s' is the stage position
                %'_t' is the time point (slidescan always has '_t1'; can it take time pts?)
                %Due to the quirks introduced in the evolution of software through ad hoc
                %programming it is some times necessary to change the format of the .tif
                %image filenames created by metamorph.
                Name_temp = regexprep(FileNames{k,j,1},'\W*',''); %Remove all not(alphabetic, numeric, or underscore) characters
                Name_temp = regexprep(Name_temp,'tocamera','','ignorecase'); %remove 'tocamera' if present b/c it is not informative
                Name_temp = regexprep(Name_temp,'camera','','ignorecase'); %remove 'camera' if present b/c it is not informative
                Name_temp = regexprep(Name_temp,'(?<=_t).*','STACK'); %'_t' is always at the end of the filename
                Name = [path2,'\Stacks\',Name_temp,'.tif']; %generates the name of the stack
                %appends the image to the tif file using no compression
                %Loading the images is the most time intensive part of the code
                name2read = [path,'\',FileNames{k,j,1}];
                newnames{j,i} = Name;
                t = Tiff(name2read,'r');
                IM = t.read;
                timelabels{i,j,k} = p53TiffMetaAnalysis4Metamorph(t);
                t.close;
                imwrite(IM,Name,'tif','WriteMode','append','Compression','none');
            end
            t = Tiff(newnames{j,i},'w');
            addTime2Stack(timelabels{i,j,:},t);
            t.close;
        end
        fprintf(1,'.'); %shows some activity to calm user... (cheap progress bar)
    end
end
fprintf(1,'\nMake Stacks Done!\n');
end

%----- SUBFUNCTION ADDTIME2STACK -----
function [] = addTime2Stack(timelabels,t)
%check for existing XML information
try
    metadata = t.getTag('ImageDescription');
    fid = fopen('t3mp.xml','w');
    fprintf(fid,'%s',metadata);
    fclose(fid);
    xdoc = xmlread('.t3mp.xml');
    node_root = xdoc.getDocumentElement;
    my_list = node_root.getElementsByTagName('time');
    aprioritimelabels = my_list.item(0).getTextContent;
    aprioritimelabels = regexp(aprioritimelabels,'.*(?=,)','match');
    timelabels = vertcat(aprioritimelabels,timelabels);
    timestr = timelabels{1};
    for i=2:length(timelabels)
        timestr = [timestr ',' timelabels{i}]; %#ok<AGROW>
    end
    my_list.item(0).setTextContent(timestr);
    xmlwrite('t3mp.xml',xdoc);
catch err %#ok<NASGU>
    timestr = timelabels{1};
    for i=2:length(timelabels)
        timestr = [timestr ',' timelabels{i}]; %#ok<AGROW>
    end
    %----- Create New XML file -----
    fid = fopen('t3mp.xml','w');
    fprintf(fid,'<MetaData></MetaData>');
    fclose(fid);
    xdoc = xmlread('.t3mp.xml');
    node_root = xdoc.getDocumentElement;
    thisElement = xdoc.createElement('time');
    thisElement.setTextContent(timestr);
    node_root.appendChild(thisElement);
    xmlwrite('t3mp.xml',xdoc);   
end
%append XML information to image description
fid = fopen('t3mp.xml');
F = fread(fid, '*char')';
fclose(fid);
t.setTag('ImageDescription',F);
end




