function [] = processManualSegTrackViaImageJ(logpath,stackpath,varargin)
% [] = processManualSegTrackViaImageJ()
% Input:
% timeReference = 'YYYYMMDD 24h:60m:60s'
%
% Output:
%
%
% Description:
%
%
% Other Notes:
% The text files that hold the manual segmentation and tracking information
% comes from ImageJ measurements. The following measurements MUST be
% included in addition to any others: "Center of Mass", "Fit Ellipse", and
% "Stack Position".
p = inputParser;
p.addRequired('logpath', @(x)ischar(x));
p.addRequired('stackpath', @(x)ischar(x));
p.addParamValue('fluorchan','YFP',@(x)ischar(x));
p.addParamValue('method','ellipse',@(x)ischar(x));
p.addParamValue('phaseratio',1,@(x)isnumeric(x));
p.parse(logpath, stackpath, varargin{:});

%% ----- Load divisions file -----
if ~isdir(logpath)
    error('manSeg:logpathNotDir','The string input for logpath is not a directory.\nPlease input the directory path, not the path to a specific file.');
end

%first, find the divisions file.
dirConLogs = dir(logpath);
temp = false;
for i=1:length(dirConLogs)
    [logFile, tok] = regexpi(dirConLogs(i).name,'divisions\.(\w+)','match','tokens','once');
    if ~isempty(logFile)
        temp = true;
        break
    end
end
if temp == false %does the division file exist
    disp(logpath)
    error('manSeg:noDiv','The %s directory does not contain a divisions file. \nPlease make sure the file is spelled with an "s".',logpath);
end
%second, check if it is MS xlsx OR just plain .txt
switch tok{1}
    %if text, read the data into a cell
    case 'txt'
        error('manSeg:divTxt','The feature to parse text files is incomplete. Please create an MS Excel file.');
        %if xlsx, collect the raw data
    case {'xls','xlsx'}
        [log_num,log_txt,~] = xlsread(fullfile(logpath,logFile));
    otherwise
        error('manSeg:divExt','The divisions file format is unrecognized. Please create a .txt or MS Excel file.');
end
%How many cell tracks are there?
numberOfCells = size(log_num,1);

%% Initialize the struct that holds all the cellular information
unitOfLife = struct('timePoints', {}, ...
    'timePoints2', {}, ...
    'time', {}, ...
    'timestamps', {}, ...
    'nucleusArea', {}, ...
    'cytoplasmArea', {}, ...
    'meanIntensity', {}, ...
    'meanIntensityHeader', {}, ...
    'parent', {}, ...
    'nuclearSolidity', {}, ...
    'divisionbool', {}, ...
    'divisionTime', {}, ...
    'manualCentroid', {}, ...
    'major', {}, ...
    'minor', {}, ...
    'angle', {}, ...
    'centroid', {}, ...
    'velocity', {}, ...
    'uid', {}, ...
    'label', {}, ...
    'originImageDirectory', {});
unitOfLife(numberOfCells).timePoints = []; %initialize the struct

%% ----- Import all of the pertinent manual segmentation and tracking
%information from a folder of text files into the unitOfLife struct. -----
pos_ind = strcmpi('position',log_txt);
pos_all = log_num(:,pos_ind);
pos_unique = unique(pos_all);
uol_ind = strcmpi('cell',log_txt);
uol_all = log_num(:,uol_ind);
dirConLogsArray = 1:length(dirConLogs);
parent_ind = strcmpi('parent',log_txt);
start_ind = strcmpi('start',log_txt);
end_ind = strcmpi('end',log_txt);
div_ind = strcmpi('division',log_txt);
div_all = num2cell(log_num(:,div_ind));
strtandend = cell(1,numberOfCells);
label_all = cell(1,numberOfCells);
for i =1:numberOfCells
    strtandend{i} = [log_num(i,start_ind) log_num(i,end_ind)];
    label_all{i} = [pos_all(i); uol_all(i)];
end
[unitOfLife(:).divisionTime] = deal(strtandend{:});
[unitOfLife(:).label] = deal(label_all{:});
parent_all = num2cell(log_num(:,parent_ind));
[unitOfLife(:).parent] = deal(parent_all{:});
[unitOfLife(:).divisionbool] = deal(div_all{:});
for j=1:numberOfCells
    %find the text file holding the segmentation and tracking info
    pos = pos_all(j);
    uol = uol_all(j);
    counter = 0;
    for k=dirConLogsArray
        counter = counter+1;
        temp = regexpi(dirConLogs(k).name,'(?<=pos)(\d+)\.(\d+)','tokens');
        if isempty(temp)
            dirConLogsArray(counter) = []; %Don't look at the same empty file more than once
            counter = counter-1;
            continue
        end
        posnum = str2double(temp{1}(1));
        uolnum = str2double(temp{1}(2));
        clear temp;
        if (posnum == pos) && (uolnum == uol)
            filename = fullfile(logpath,dirConLogs(k).name);
            manualData = importdata(filename);
            dataOffset = size(manualData.textdata,2)-size(manualData.data,2);
            headers={'XM';'YM';'Major';'Minor';'Angle';'Slice'}; %BEWARE! BEFORE CHANGING THE ORDER OF THESE HEADERS RECOGNIZE THE CONSEQUENCES.
            index=zeros(size(headers)); %The column numbers for the wanted information will be stored here
            for i=1:length(headers)
                temp = strcmp(headers(i),manualData.textdata(1,:));
                if any(temp)
                    index(i) = find(temp,1,'first'); %get labels
                elseif strcmp(headers{i},'Slice') %if there are no slices the same information can be extracted from the label should it exist
                    temp = strcmp('Label',manualData.textdata(1,:));
                    if any(temp)
                        index(i) = find(temp,1,'first'); %get labels
                        index2 = size(manualData.data,2)+1;
                        for h = 2:size(manualData.textdata,1)
                            temp = regexp(manualData.textdata{h,index(i)},'(?<=:)\d+','match');
                            manualData.data(h-1,index2) = str2double(temp);
                        end
                        index(i)=index2+dataOffset;
                    else
                        error('manSeg:missingData2','The "Label" and "Slice" data are missing from the manual segmentation and tracking text file.');
                    end
                else %there is data missing
                    error('manSeg:missingData1','The "%s" data is missing from the manual segmentation and tracking text file.',headers{i});
                end
            end
            %As an idiosyncrasy of using importdata() with the ImageJ
            %text files the data must be parsed in a specific manner.
            %Therefore, there is code written to do this parsing that may
            %also be a hotspot for bugs.
            index = index - dataOffset; %The numeric data will not contain columns with text data and all of the text data is found in the first column(s). All columns have text headers, so the index of the text column(s) needs to be removed.
            %Fill the unitOfLife struct with the relevant data from the
            %text file
            unitOfLife(j).manualCentroid = [manualData.data(:,index(2)),manualData.data(:,index(1))];
            unitOfLife(j).major = manualData.data(:,index(3));
            unitOfLife(j).minor = manualData.data(:,index(4));
            unitOfLife(j).angle = manualData.data(:,index(5));
            unitOfLife(j).timePoints = manualData.data(:,index(6));
            %adjust the timepoints, derived from the phase image, to match the number of timepoints found in the fluorescent images. For example, the phase images might be capture every 5 minutes, whereas the fluorescent images every 20 minutes.
            timePointsTemp = (unitOfLife(j).timePoints + p.Results.phaseratio - 1)/p.Results.phaseratio;
            timePointsLogical = mod(timePointsTemp,1)==0;
            timePointsTemp = timePointsTemp(timePointsLogical);
            unitOfLife(j).timePoints = timePointsTemp;
            unitOfLife(j).manualCentroid = unitOfLife(j).manualCentroid(timePointsLogical,:);
            unitOfLife(j).major = unitOfLife(j).major(timePointsLogical);
            unitOfLife(j).minor = unitOfLife(j).minor(timePointsLogical);
            unitOfLife(j).angle = unitOfLife(j).angle(timePointsLogical);
            dirConLogsArray(counter) = []; %Don't look at the same text file more than once
            break
        end
    end
end

%% ----- Import PNG paths -----
%Check whether the computer is a mac
if ismac
    error('manSeg:notApc','Flatfield correction is currently configured to run on a pc');
    %[~,filepaths]=system('ls -Rp */*.png');
elseif ispc
    %the 'dir' command in Matlab does not search subfolders. However,
    %there is a MS-DOS command that will do the trick!
    %'/B' uses bare format (no heading information or summary)
    %'/S' displays files in specified directory and all subdirectories
    cd(stackpath)
    [~,filepaths] = system('dir /S /B *.png');
    filepaths = textscan(filepaths,'%s');
    filepaths = filepaths{1};
else
    error('manSeg:notApc','processManualSegTrackViaImageJ.m is currently configured to run on a pc');
end

%% Find numbers needed for data extraction from images and create a "logic-space" that will guide the data extraction
%Find the ultimate timepoint that was manually segmented
numberOfTimepoints = 0;
for i=1:numberOfCells
    timemax4acell = max(unitOfLife(i).timePoints);
    if timemax4acell>numberOfTimepoints
        numberOfTimepoints = timemax4acell;
    end
end
numberOfPositions = length(pos_unique);
%With the number of positions, the number of cells in each position, and
%the number of timepoints, create a 3 dimensional logic-space that will
%contain a 1 if there is a cell in a particular position at a particular
%time and a 0 otherwise. This logic space will guide the operation of
%opening an image and extracting the data for all the cells contained
%within that image.
logicspace = zeros(numberOfCells, numberOfTimepoints, numberOfPositions);
for i=1:numberOfCells
    logicspace(i,unitOfLife(i).timePoints,pos_unique == pos_all(i)) = true;
end

%% Extract the data from the images
for i=1:numberOfPositions
    %find the subset of paths which are images from a particular position
    filepaths4posi = cell(numberOfTimepoints,1);
    for j=1:length(filepaths)
        tok = regexpi(filepaths{j},'_s(\d+)_w\d+(.+)_t(\d+)','tokens');
        if (~isempty(tok)) && (strcmp(tok{1}{2},p.Results.fluorchan)) && (str2double(tok{1}{1}) == pos_unique(i))
            filepaths4posi{str2double(tok{1}{3})} = filepaths{j};
        end
    end
    for j=1:numberOfTimepoints
        %extract data using the ellipse method
        if ~isempty(filepaths4posi{j})
            IM = imread(filepaths4posi{j});
            [unitOfLife]=distillDataFromMovie(logicspace(:,j,i),unitOfLife,'ellipse',IM,j);
        end
    end
end
filename = sprintf('dynamics%s',p.Results.fluorchan);
save(fullfile(logpath,filename),'unitOfLife');
%% link the cell segments together
linkmap = cell(1,length(unitOfLife));
parent_all = [unitOfLife(:).parent];
poslabel_all = [unitOfLife(:).label];
div_all = [unitOfLife(:).divisionbool];
rootnodepointer = 1;
linkcounter = 1;
currentbranch = [];
childremovedflag = false;
%find all the root nodes
rootnodes = parent_all == 0;
rootnodes = find(rootnodes);
realrootnodeslength = length(rootnodes);
rootnodes(end+1) = 0;
cellpointer = rootnodes(rootnodepointer);
while rootnodepointer <= realrootnodeslength
    currentbranch(end + 1) = cellpointer; %#ok<AGROW>
    if div_all(cellpointer) == true
        %Where are the child nodes
        posid = poslabel_all(1,cellpointer);
        cellid = poslabel_all(2,cellpointer);
        temp1 = poslabel_all(1,:)==posid;
        temp2 = parent_all==cellid;
        temp3 = temp1 & temp2;
        if any(temp3)
            cellpointer = find(temp3,1,'first');
        else
            if childremovedflag
                if any(rootnodes == cellpointer)
                    rootnodepointer = rootnodepointer + 1;
                    cellpointer = rootnodes(rootnodepointer);
                    currentbranch = [];
                    childremovedflag = false;
                else
                    parent_all(cellpointer) = 0;
                    currentbranch(end) = [];
                    cellpointer = currentbranch(end);
                    currentbranch(end) = [];
                end
            else
                if any(rootnodes == cellpointer)
                    linkmap{linkcounter} = currentbranch;
                    linkcounter = linkcounter + 1;
                    rootnodepointer = rootnodepointer + 1;
                    cellpointer = rootnodes(rootnodepointer);
                    currentbranch = [];
                else
                    linkmap{linkcounter} = currentbranch;
                    linkcounter = linkcounter + 1;
                    parent_all(cellpointer) = 0;
                    currentbranch(end) = [];
                    cellpointer = currentbranch(end);
                    currentbranch(end) = [];
                    childremovedflag = true;
                end
            end
        end
    else
        if any(rootnodes == cellpointer)
            linkmap{linkcounter} = currentbranch;
            linkcounter = linkcounter + 1;
            rootnodepointer = rootnodepointer + 1;
            cellpointer = rootnodes(rootnodepointer);
            currentbranch = [];
        else
            linkmap{linkcounter} = currentbranch;
            linkcounter = linkcounter + 1;
            parent_all(cellpointer) = 0;
            currentbranch(end) = [];
            cellpointer = currentbranch(end);
            currentbranch(end) = [];
            childremovedflag = true;
        end
    end
end
linkmap(linkcounter:end) = [];
disp('cool')

function [unitOfLife]=distillDataFromMovie(map,unitOfLife,method,IM,tind)
%% The methods available for gathering single cell measurements are:
%1. ellipse: The shape used for segmentation is an ellipse.
switch lower(method)
    case 'ellipse'
        %Initialize the arrays that will contain the points that describe
        %the ellipse
        xellipse=zeros(38,1);
        yellipse=zeros(38,1);
        rho = (0:9:333)/53;
        rhocos=cos(rho);
        rhosin=sin(rho);
        mapIndex = find(map)';
        for i=mapIndex
            k = unitOfLife(i).timePoints == tind;
            if sum(k) > 1
                msg = sprintf('\n%d at time %d has a problem with the log file entries\n',i, tind);
                disp(msg);
                continue
            end
            %All numbers that describe the ellipse must be rounded in
            %order to be discretized, i.e. refer to a pixel coordinate
            if unitOfLife(i).angle(k) == 0 || unitOfLife(i).angle(k) == 180
                a = round(unitOfLife(i).major(k)/2);
                b = round(unitOfLife(i).minor(k)/2);
            elseif unitOfLife(i).angle(k) == 90 || unitOfLife(i).angle(k) == 270
                a = round(unitOfLife(i).minor(k)/2);
                b = round(unitOfLife(i).major(k)/2);
            end
            %Centroid is determined
            xm = round(unitOfLife(i).manualCentroid(k,1));
            ym = round(unitOfLife(i).manualCentroid(k,2));
            
            %the points on the perimeter of the ellipse with the determined
            %centroid and major/minor axes are calculated
            for k=1:38
                x=xm+a*rhocos(k);
                y=ym+b*rhosin(k);
                xellipse(k) = round(x);
                yellipse(k) = round(y);
            end
            %The pixels within the elipse are determined with the function
            %roipoly and their mean intensities are stored in the output array
            BW = roipoly(IM,yellipse,xellipse); %BW is a binary mask
            IM_temp = IM;
            unitOfLife(i).meanIntensity(end+1) = mean(IM_temp(BW));
            unitOfLife(i).timePoints2(end+1) = tind;
            fprintf(1,'.');
        end
    case 'watershed'
    otherwise
        error('unknown segmentation method input into getMovieMeasurements.m')
end
