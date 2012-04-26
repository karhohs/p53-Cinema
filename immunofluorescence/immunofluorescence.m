function [] = immunofluorescenceCP(IFpath,varargin)
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
% This program will segment nuclei and cytoplasm if markers for these
% exist. It is assumed that at least a nuclear marker exists.
p = inputParser;
p.addRequired('IFpath', @(x)ischar(x));
p.addParamValue('nuclearMarker','DAPI',@(x)ischar(x));
p.addParamValue('timeReference','',@(x)ischar(x)||isnumeric(x));
p.addParamValue('timeUnits','hours',@(x)ischar(x));
p.addParamValue('method','ellipse',@(x)ischar(x));
p.parse(logpath, stackpath, varargin{:});

%Store all of the biological information in a struct
unitOfLife = struct('timePoints', {}, 'time', {}, 'timeUnits', {}, 'nucleusArea', {}, 'cytoplasmArea', {}, 'meanIntensity', {},'meanIntensityHeader', {},'parent', {}, 'nuclearSolidity', {}, 'divisionTime', {}, 'manualCentroid', {}, 'major', {},'minor', {}, 'angle', {},'centroid', {},'velocity', {}, 'uid', {}, 'originImage', {});
unitOfLife(1e8).time = []; %initialize the struct
cellcounter = 1; %keeps track of the number of cells counted


%Load the directory contents with the image files. Look for the DAPI
%channel.
dirCon_stack = dir(stackpath);
    %Find the stack of images of the specified channel
    stack_filename = '';
    for j=1:length(dirCon_stack)
        expr = '.*_w\d+(.+)_s(\d+).*';
        temp = regexp(dirCon_stack(j).name,expr,'tokens');
%identify the channels
    expr='(?<=_w\d+).*(?=_s\d+)';

%find the time of the image
t = Tiff(name2read,'r');
timelabels{k} = p53TiffMetaAnalysis4Metamorph(t);
t.close;
                
end