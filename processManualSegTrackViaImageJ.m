function [reltime,time,meanGreyVal] = processManualSegTrackViaImageJ(logpath,stackpath,varargin)
% [] = processManualSegTrackViaImageJ()
% Input:
%
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
p.addParamValue('timeReference','',@(x)ischar(x));
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
unitOfLife = struct('timePoints', {}, 'time', {}, 'nucleusArea', {}, 'cytoplasmArea', {}, 'meanIntensity', {},'parent', {}, 'nuclearSolidity', {}, 'divisionTime', {}, 'manualCentroid', {}, 'major', {},'minor', {}, 'angle', {},'centroid', {},'velocity', {}, 'uid', {}, 'originImage', {});
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
%It is convenient to have a matrix that helps identify which cells are
%present at each timepoint.
%The stacks to open are determined by the unique positions
%for each position...
for i=1:length(pos_unique)
    %Find the stack of images of the specified channel
    dirCon_stack = dir(stackpath);
    for j=1:length(dirCon_stack)
        expr = '.*_w\d+(.+)_s(\d+).*';
        temp = regexp(dirCon_stack(j).name,expr,'tokens');
            if (~isempty(temp)) && (strcmp(temp{1}{1},p.Results.fluorchan)) && (str2double(temp{1}{2}) == pos_unique(i))
                stack_filename = dirCon_stack(j).name;
                break
            end
    end
    %Import stack of images
    filename_stack = fullfile(stackpath,stack_filename);
    t = Tiff(filename_stack,'r');
    info = imfinfo(filename_stack,'tif');
    %Import the time information
    timeTable = importTimeTableFromStack(t);
    %Create a time vector (units = hours) relative to the timeref
    relTimeTable = createRelTimeTable(timeTable);
    %Create a data structure that will hold the tracking info for all the cells in the current position.
    trackingTable = cell(length(ceTable{i}),length(info));
    %Import the tracking data for each cell in this stage position
    for j=1:length(ceTable{i})
        %find the text file holding the segmentation and tracking info
        for k=1:length(dirConLogs)
            temp = regexp(dirConLogs(k).name,'(?<=Pos)(\d+)\.(\d+)','tokens');
            pnum = str2num(temp{1});
            cenum = str2num(temp{2});
            if pnum == position(i) && cenum == ceTable{i}(j)
                filename_log = fullfile(logpath,dirConLogs(k).name);
                break
            end
        end
        %import this data
        trackingTable(j,:) = importTrackingData(filename_log);
    end
    %for each cell in that stack...
    for j=1:length(info)
        for k=1:length(ceTable{i})
            %distill the data from the movie
            distillDataFromMovie();
            %store this data in a MATLAB-cell and a struct, two different ways to store the data
            
        end
    end
    t.close;
end
%Create a struct as an alternative representatin of the data
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