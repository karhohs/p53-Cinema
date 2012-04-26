function [] = immunofluorescenceCP(CPoutput,varargin)
% [] = processManualSegTrackViaImageJ()
% Input:
% CPoutput = the path to a Cell Profiler *.mat file.
%
% Output:
%
%
% Description:
% I have been liking Cell Profiler for quantifying immunofluorescence data
% recently. This file will analyze the output from Cell Profiler.
%
% Other Notes:
% 
p = inputParser;
p.addRequired('CPoutput', @(x)ischar(x));
p.addParamValue('infoName','Nuclei',@(x)ischar(x));
p.parse(CPoutput, varargin{:});

Nuclei = p.Results.infoName;
CP = open(CPoutput);

%Plot some histograms
nucleiArea = eval(['CP.handles.Measurements.' Nuclei '.AreaShape_Area']);
nucleiMeanIntensityYFP = eval(['CP.handles.Measurements.' Nuclei '.Intensity_MeanIntensity_YFP']);
nucleiMeanIntensityCy5 = eval(['CP.handles.Measurements.' Nuclei '.Intensity_MeanIntensity_Cy5']);
nucleiMeanIntensityDAPI = eval(['CP.handles.Measurements.' Nuclei '.Intensity_MeanIntensity_DAPI']);
nucleiSolidity = eval(['CP.handles.Measurements.' Nuclei '.AreaShape_Solidity']);
nucleiCompactness = eval(['CP.handles.Measurements.' Nuclei '.AreaShape_Compactness']);
nucleiX = eval(['CP.handles.Measurements.' Nuclei '.AreaShape_Center_Y']);
nucleiY = eval(['CP.handles.Measurements.' Nuclei '.AreaShape_Center_X']);
DAPIfilenames = CP.handles.Measurements.Image.FileName_DAPI;

numberoflifeunits = 0;
for i = 1:length(nucleiArea)
    numberoflifeunits = numberoflifeunits + numel(nucleiArea{i});
end

nucleiArea_array = linearizeContentsOfCell(nucleiArea,numberoflifeunits);
nucleiMeanIntensityYFP_array = linearizeContentsOfCell(nucleiMeanIntensityYFP,numberoflifeunits);
nucleiMeanIntensityCy5_array = linearizeContentsOfCell(nucleiMeanIntensityCy5,numberoflifeunits);
nucleiMeanIntensityDAPI_array = linearizeContentsOfCell(nucleiMeanIntensityDAPI,numberoflifeunits);
nucleiSolidity_array = linearizeContentsOfCell(nucleiSolidity,numberoflifeunits);
nucleiCompactness_array = linearizeContentsOfCell(nucleiCompactness,numberoflifeunits);

%----- Create a new folder to hold corrected images -----
tempfoldername=regexp(CPoutput,'(?<=\\)[\w ]*','match'); %Prepare to create a new folder to place background subtracted stacks
tempdrive = regexp(CPoutput,'[a-zA-Z]:','match');
tempfoldername{end} = 'figures';
imagepath=fullfile(tempdrive{1},tempfoldername{:});
mkdir(imagepath);
cd(imagepath)

myplothist(nucleiArea_array,'nucleiArea')
myplothist(nucleiSolidity_array,'nucleiSolidity')
myplothist(nucleiCompactness_array,'nucleiCompactness')
myplothist(nucleiMeanIntensityDAPI_array.*nucleiArea_array,'nucleiDAPI')
end

function [y] = linearizeContentsOfCell(x,numolu)
y = zeros(1,numolu);
counter = 1;
for i = 1:length(x)
    temp = numel(x{i});
    y(counter:temp+counter-1) = x{i};
    counter = counter+temp;
end
end

function [] = myplothist(x,filename)
h=dialog ( 'visible', 'off', 'windowstyle', 'normal' );
ax=axes('parent', h);
hist (ax,x,100)
title(filename)
%plot (ax,[1 2 3 4 5], [1 2 3 4 5].^3 )
saveas (ax, filename, 'jpg' )
close(h)
end