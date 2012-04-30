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

CP = open(CPoutput);

%Plot some histograms
nucleiArea = CP.handles.Measurements.Nuclei.AreaShape_Area;
nucleiMeanIntensityYFP = CP.handles.Measurements.Nuclei.Intensity_MeanIntensity_YFP;
nucleiMeanIntensityCy5 = CP.handles.Measurements.Nuclei.Intensity_MeanIntensity_Cy5;
nucleiMeanIntensityDAPI = CP.handles.Measurements.Nuclei.Intensity_MeanIntensity_DAPI;
nucleiSolidity = CP.handles.Measurements.Nuclei.AreaShape_Solidity;
nucleiCompactness = CP.handles.Measurements.Nuclei.AreaShape_Compactness;
nucleiFormFactor = CP.handles.Measurements.Nuclei.AreaShape_FormFactor;
nucleiX = CP.handles.Measurements.Nuclei.AreaShape_Center_Y;
nucleiY = CP.handles.Measurements.Nuclei.AreaShape_Center_X;
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
nucleiFormFactor_array = linearizeContentsOfCell(nucleiFormFactor,numberoflifeunits);

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
myplothist(nucleiFormFactor_array,'nucleiFormFactor')
myplothist(nucleiMeanIntensityDAPI_array.*nucleiArea_array,'TotalNucleiDAPI')
myplothist(nucleiMeanIntensityYFP_array.*nucleiArea_array,'TotalNucleiYFP')
myplothist(nucleiMeanIntensityCy5_array.*nucleiArea_array,'TotalNucleiCy5')

%Filtering the dataset
%Area is used to filter out objects that are too small to be nuclei. Use a
%log-normal distribution to eliminate small outliers. Senescent cells are
%known to have large nuclei, so these outliers should not be eliminated.
afarray = areafilter()


%Clumped cells and senescent cells have similar areas. To distinguish
%between clumped cells and senescent cells a combined measure is used:
%compactness, form-factor, and solidity. These 3 measures give a sense of
%how circular, or round, is a nuclei. Senescent cells tend to have very
%large round nuclei.
%form-factor gives a measure of how circular an object is by comparing the
%equivalent area of circle to the perimenter. This measure 
%should be less than 1. It is best fit to a distribution when the measures
%are abs(1-x), where x is the form-factor measurement. Since 0<x<1 it might
%be worth looking at the -log(x) distribution.

%Solidity detects the presence of furrows or deviations from a smooth curve
%by comparing the convex hull area to the object area. Clumped cells can be
%thought of a overlapping ellipses that would produce object furrows. It is
%a measure that is 0<x<1 and has a similar distribution to the form-factor.

%Compactness gives a second way to measure furrows or spikes.
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
h=figure ( 'visible', 'off', 'position', [10, 10, 672, 512] );
ax=axes('parent',h);
hist(ax,x,100);
title(filename);
saveas (ax, filename, 'jpg' );
close(h);
end

function [y] = areafilter(x)
L = length(x);
%assume area is distributed as a log-normal distribution
logx = log(x);
%assume there are outliers. We want to ignore these outliers while fitting
%the normal distribution, so we discard the bottom and top 10%.
logx = logx(round(L*.1):round(L*.9));
%From this trimmed set we have to bootstrap the tails of the distribution.
bootset = randn(L,1);
bootset = sort(bootset)
%We must estimate the variance of the distribution from the trimmed set.
%This variance must be scaled by the variance of a zero mean, unit variance
%normal distribution trimmed the same way. The z-score for a 10% tail is
%1.2815.
%The variance for this trimmed distribution can be found numerically using
%MATLAB. x = randn(1e6,1);x = sort(x);x = x(round(length(x)*.1):round(length(x)*.9));var(x)
%he variance equals 0.4377
varlogx = var(logx)/.4377;
%The bootset is scaled by the mean and variance of the sample data
bootset = logx*sqrt(varlogx)+mean(logx);
%The bootstraped tails are added to the trimmed distribution.
bootmin = bootset(bootset<min(logx));
bootmax = bootset(bootset>max(logx));
newdist = [bootmin;logx';bootmax];
%The new distribution is scaled to zero mean and unit variance
newdistmean = mean(newdist);
newdistvar = var(newdist);
newdist = (newdist-newdistmean)/sqrt(newdistvar);
%The goodness of fit to a normal distribution is tested by the
%kolmogorov-smirnoff test.
goodfitbool = kstest(newdist);
end