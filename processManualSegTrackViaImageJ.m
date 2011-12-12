function [reltime,time,meanGreyVal] = processManualSegTrackViaImageJ(logpath,stackpath,channel,timeref)
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

%Load divisions file
%first, find the divisions file.
dirConLogs = dir(logpath);
for i=1:length(dirConLogs)
    [mat, tok] = regexp(dirConLogs(i).name,'divisions\.(\w+)','match','tokens');
    if ~isempty(mat)
        temp = true;
        break
    end
end
if temp == false
    error('manSeg:noDiv',['The ', logpath, ' directory does not contain a divisions file']);
end
%second, check if it is MS xlsx OR just plain .txt
switch tok
    %if text, read the data into a cell
    case 'txt'
        error('manSeg:divTxt','The feature to parse text files is incomplete. Please create a MS Excel file.');
        %if xlsx, collect the raw data
    case {'xls','xlsx'}
        [div_num,div_txt,~] = xlsread();
    otherwise
        error('manSeg:divExt','The divisions file format is unrecognized. Please create a .txt or MS Excel file.');
end
%Identify the stacks of images that contain the cells of interest
pos_ind = strcmpi('position',div_txt);
position_temp = div_num(:,pos_ind);
position = unique(position_temp);
%Sort the cells by their image origin into a MATLAB-cell
ceTable = cell(size(position));
ce_ind = strcmpi('cell',div_txt);
ce_temp = div_num(:,ce_ind);
for i=1:length(ceTable)
    ceTable{i} = ce_temp(position_temp==position(i));
end
%initialize outputs
time = cell(size(position));
reltime = cell(size(position));
meanGreyVal =cell(size(position));
%for each stack...
for i=1:length(position)
    %Find the stack of images of the specified channel
    dirCon_stack = dir(stackpath);
    position_str = num2str(position(i));
    for j=1:length(dirCon_stack)
        temp = regexp(dirCon_stack(i).name,channel,'once');
        if ~isempty(temp)
            temp = regexp(dirCon_stack(i).name,'_s(\d+)','tokens');
            if temp == position_str
                stack_filename = dirCon_stack(i).name;
                break
            end
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