function [reltime,time,meanGreyVal] = processManualSegTrackViaImageJ(logpath,stackpath,varargin)
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
p.addParamValue('timeReference','',@(x)ischar(x)||isnumeric(x));
p.addParamValue('timeUnits','hours',@(x)ischar(x));
p.parse(logpath, stackpath, varargin{:});
%----- Load divisions file -----
%first, find the divisions file.
dirConLogs = dir(logpath);
for i=1:length(dirConLogs)
    [logFile, tok] = regexpi(dirConLogs(i).name,'divisions\.(\w+)','match','tokens','once');
    if ~isempty(logFile)
        temp = true;
        break
    end
end
if temp == false %does the division file exist
    error('manSeg:noDiv',['The ', logpath, ' directory does not contain a divisions file']);
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
%Initialize the struct that holds all the cellular information
unitOfLife = struct('timePoints', {}, 'time', {}, 'timeUnits', {}, 'nucleusArea', {}, 'cytoplasmArea', {}, 'meanIntensity', {},'parent', {}, 'nuclearSolidity', {}, 'divisionTime', {}, 'manualCentroid', {}, 'major', {},'minor', {}, 'angle', {},'centroid', {},'velocity', {}, 'uid', {}, 'originImage', {});
unitOfLife(numberOfCells).time = []; %initialize the struct
%----- Import all of the pertinent manual segmentation and tracking
%information from a folder of text files into the unitOfLife struct. -----
pos_ind = strcmpi('position',log_txt);
pos_all = log_num(:,pos_ind);
pos_unique = unique(pos_all);
uol_ind = strcmpi('cell',log_txt);
uol_all = log_num(:,uol_ind);
dirConLogsArray = 1:length(dirConLogs);
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
            headers={'XM';'YM';'Major';'Minor';'Angle';'Slice'}; %BEWARE! BEFORE CHANGING THE ORDER OF THESE HEADERS RECOGNIZE THE CONSEQUENCES.
            index=zeros(size(headers)); %The column numbers for the wanted information will be stored here
            for i=1:length(headers)
                index(i) = find(strcmp(headers(i),manualData.textdata(1,:)),1,'first'); %get labels
            end
            %As an idiosyncrasy of using importdata() with the ImageJ
            %text files the data must be parsed in a specific manner.
            %Therefore, there is code written to do this parsing that may
            %also be a hotspot for bugs.
            dataOffset = size(manualData.textdata,2)-size(manualData.data,2);
            index = index - dataOffset; %The numeric data will not contain columns with text data and all of the text data is found in the first column(s). All columns have text headers, so the index of the text column(s) needs to be removed.
            %Fill the unitOfLife struct with the relevant data from the
            %text file
            unitOfLife(j).manualCentroid = [manualData.data(:,index(2)),manualData.data(:,index(1))];
            unitOfLife(j).major = manualData.data(:,index(3));
            unitOfLife(j).minor = manualData.data(:,index(4));
            unitOfLife(j).angle = manualData.data(:,index(5));
            unitOfLife(j).timePoints = manualData.data(:,index(6));
            dirConLogsArray(counter) = []; %Don't look at the same text file more than once
            break
        end
    end
end
%Our goal is to only open each stack of images once.
%The stacks to open are determined by the unique positions
%for each position...
dirCon_stack = dir(stackpath);
for i=1:length(pos_unique)
    %Find the stack of images of the specified channel
    stack_filename = '';
    for j=1:length(dirCon_stack)
        expr = '.*_w\d+(.+)_s(\d+).*';
        temp = regexp(dirCon_stack(j).name,expr,'tokens');
        if (~isempty(temp)) && (strcmp(temp{1}{1},p.Results.fluorchan)) && (str2double(temp{1}{2}) == pos_unique(i))
            stack_filename = dirCon_stack(j).name;
            break
        end
    end
    if isempty(stack_filename)
        warning('ManSeg:stackMissing', ...
            'Could not locate stack for position "%s".', pos_unique(i))
        continue
    end
    %Import stack of images (this is probably done inefficiently)
    filename_stack = fullfile(stackpath,stack_filename);
    info = imfinfo(filename_stack,'tif');
    sizeOfImage = [info(1).Height, info(1).Width, length(info)];
    IM = zeros(sizeOfImage);
    IM = loadStack(filename_stack,IM,sizeOfImage);
    %Import the time information
    [time,tUnit] = importTimeFromStack(filename_stack,p.Results.timeReference,p.Results.timeUnits);
    %It is convenient to have a matrix that helps identify which cells are
    %present at each timepoint.
    numberOfCellsArray = 1:numberOfCells;
    numberOfCellsArray = numberOfCellsArray(pos_all==pos_unique(i));
    map = zeros(length(numberOfCellsArray),length(time{3}));
    %create the matrix and a unique ID for each cell and add the time and
    %units
    for k=1:length(numberOfCellsArray)
        j = numberOfCellsArray(k);
        unitOfLife(j).time = time{3}(unitOfLife(j).timePoints);
        unitOfLife(j).timeUnits = tUnit;
        map(k,unitOfLife(j).timePoints) = j;
        unitOfLife(j).uid = ['pos ' num2str(pos_all(j)) ' cell ' num2str(uol_all(j)) ' timestamp ' time{1}{unitOfLife(j).timePoints(1)}];
        unitOfLife(j).originImage = [filename_stack ':' num2str(unitOfLife(j).timePoints(1))];
    end
    %for each cell in that stack...
    for j=1
        for k=1
            %distill the data from the movie
            distillDataFromMovie();
            
            
        end
    end
end
end

function [data]=distillDataFromMovie(segment,method,channels,path)
%% The methods available for gathering single cell measurements are:
%1. ellipse: The shape used for segmentation is an ellipse.

switch lower(method)
    case 'ellipse'
        %Find the max time
        Tmax=0;
        for i=1:length(segment)
            a=max(segment(i).time);
            if a>Tmax
                Tmax=a;
            end
        end
        data=struct;
        %Initialize the arrays that will contain the points that describe
        %the ellipse
        xellipse=zeros(38,1);
        yellipse=zeros(38,1);
        rho = (0:9:333)/53;
        rhocos=cos(rho);
        rhosin=sin(rho);
        
        % Remove 'Brightfield Camera', since no protein measurements are
        % found in this channel
        temp=cell(1,length(channels)-1);
        i=1;
        for h=1:length(channels)
            if isempty(regexp(channels{h},'Phase|Bright\w*','once'))
                temp{i}=channels{h};
                i=i+1;
                %else
                %bf=h;
            end
        end
        channels2=temp;
        clear temp
        
        for h=1:length(channels2)
            
            for i=1:length(segment)
                fprintf('\n'); disp(['crunching data in ',channels2{h},' of Pos',segment(i).pos])
                %Load the stack for a given position and fluorscent wavelength
                fname=[path,'\Pos',segment(i).pos(1:2),channels2{h},'.TIF'];
                S=readStack(fname);
                temp=NaN(Tmax,1);
                for j=1:length(segment(i).time)
                    %All numbers that describe the ellipse must be rounded in
                    %order to be discretized, i.e. refer to a pixel coordinate
                    if segment(i).ang(j) == 0 || segment(i).ang(j) == 180
                        a = round(segment(i).maj(j)/2);
                        b = round(segment(i).mnr(j)/2);
                    elseif segment(i).ang(j) == 90 || segment(i).ang(j) == 270
                        a = round(segment(i).mnr(j)/2);
                        b = round(segment(i).maj(j)/2);
                    end
                    %Centroid is determined
                    xm = round(segment(i).xm(j));
                    ym = round(segment(i).ym(j));
                    
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
                    BW = roipoly(S{segment(i).time(j)},xellipse,yellipse); %BW is a binary mask
                    temp(segment(i).time(j)) = mean(mean(S{segment(i).time(j)}(BW)));
                    fprintf(1,'.');
                end
                %Remove the spaces from channel name
                %This allows the channel names to be used as field names
                name=regexprep(channels2{h},' ','');
                data(i).(name)=temp;
                data(i).pos=segment(i).pos;
            end
        end
        
    case 'watershed'
    otherwise
        error('unknown segmentation method input into getMovieMeasurements.m')
end
end

function [IM] = loadStack(path,IM,s)
t = Tiff(path,'r');
if s(3) > 1
    for k=1:s(3)-1
        IM(:,:,k) = double(t.read);
        t.nextDirectory;
    end
end
%one last time without t.nextDirectory
IM (:,:,s(3)) = double(t.read);
t.close;
end

function [time,tUnits] = importTimeFromStack(filename,tRef,tUnits)
t = Tiff(filename,'r');
%import the time metadata
metadata = t.getTag('ImageDescription');
fid = fopen('t3mp.xml','w');
fprintf(fid,'%s',metadata);
fclose(fid);
xdoc = my_parseXML('t3mp.xml');
delete('t3mp.xml');
%import times into a cell of strings. The text strings are comma separated
%values.
time = textscan(xdoc.MetaData.time.txtd8a{1},'%s','delimiter', ',');
numtimepts = length(time{1});
time{2} = 1:numtimepts;
%make sure tRef is a valid input
if isempty(tRef)
    tRef = time{1}{1};
elseif ischar(tRef)
    
elseif isnumeric(tRef)
    if (rem(tRef,1)~=0) || (tRef<1) || (tRef>numtimepts)
        warning('ManSeg:tRef', ...
            'tRef, %s, was an invalid time point and was ignored', tRef)
        tRef = time{1}{1};
    else
        tRef = time{1}{tRef};
    end
end
%convert time measurements in relative time measurements at the tUnits
%scale. The time is saved in the following format: 'YYYYMMDD 24h:60m:60s'
%repeat the loop once more for the time reference
C = regexp(tRef,'(\d{4})(\d{2})(\d{2}) (\d{2}):(\d{2}):(.+)','tokens');
tRef = zeros(1,6);
if ~isempty(C)
    tRef(1) = str2double(C{1}{1});
    tRef(2) = str2double(C{1}{2});
    tRef(3) = str2double(C{1}{3});
    tRef(4) = str2double(C{1}{4});
    tRef(5) = str2double(C{1}{5});
    tRef(6) = str2double(C{1}{6});
else
    error('manSeg:timeFormat','The expected time format was not found.');
end
%calculate the difference in time
tRef = datenum(tRef);
DateVector = zeros(1,6);
for i=1:numtimepts
    C = regexp(time{1}{i},'(\d{4})(\d{2})(\d{2}) (\d{2}):(\d{2}):(.+)','tokens');
    if ~isempty(C)
        DateVector(1) = str2double(C{1}{1});
        DateVector(2) = str2double(C{1}{2});
        DateVector(3) = str2double(C{1}{3});
        DateVector(4) = str2double(C{1}{4});
        DateVector(5) = str2double(C{1}{5});
        DateVector(6) = str2double(C{1}{6});
        tDiff = datenum(DateVector) - tRef; %This number is in units of days. The key number is 86400 seconds in a day
        switch lower(tUnits)
            case {'second','seconds'}
                time{3}(i) = tDiff*86400;
                tUnits = 'seconds';
            case {'minute','minutes'}
                time{3}(i) = tDiff*1440;
                tUnits = 'minutes';
            case {'hour','hours'}
                time{3}(i) = tDiff*24;
                tUnits = 'hours';
            case {'day','days'}
                time{3}(i) = tDiff;
                tUnits = 'days';
            case {'week','weeks'}
                time{3}(i) = tDiff/7;
                tUnits = 'weeks';
            case {'month','months'}
                time{3}(i) = tDiff/30.4167;
                tUnits = 'months';
            case {'year','years'}
                time{3}(i) = tDiff/365;
                tUnits = 'years';
            otherwise
                warning('ManSeg:tUnits', ...
                    'tUnits, %s, was not recognized and converted into "hours"', tUnits)
                time{3}(i) = tDiff*24;
                tUnits = 'hours';
        end
    else
        error('manSeg:timeFormat','The expected time format was not found.');
    end
end
end