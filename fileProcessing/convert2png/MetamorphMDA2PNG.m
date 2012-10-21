function [] = MetamorphMDA2PNG(path)
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
path2=pwd;
warning('off','MATLAB:tifflib:libraryWarning');
% ----- Organizing the filenames -----
%The filenames are organized within a 3D (i,j,k) matrix where each dimension
%represents either time, position, or wavelength (i=time, j=position,
%k=wavelength). This organization facilitates the collating these images
%into TIFF stacks.

FileNames=cell(2048,256,8); %Assumes there will not be more than 2048 timepoints (over a week at images every 5 minutes), 256 positions, or more than 8 channels. If there are, change this number. It has been pre-allocated for speed.
%FileNames(Time,Position,Channel)
disp(['working in ', path]); %Collect the directory contents
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
Tunique = unique(T);
Tunique = Tunique(Tunique~=0);
Sunique = unique(S);
Sunique = Sunique(Sunique~=0);
Wunique = unique(W);
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
%create the root directory that will store the PNG images

disp('cool');