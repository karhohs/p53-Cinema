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


%Filtering the dataset
%Area is used to filter out objects that are too small to be nuclei. Use a
%log-normal distribution to eliminate small outliers. Senescent cells are
%known to have large nuclei, so these outliers should not be eliminated.
[afarray,afthresh,mu,sigma] = areafilter(nucleiArea_array);
myplothistarea(nucleiArea_array,'nucleiArea',afthresh,mu,sigma);

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
[fffarray,fffthresh] = formfactorfilter(nucleiFormFactor_array);
myplothistformfactor(nucleiFormFactor_array,'nucleiFormFactor',fffthresh)
%Solidity detects the presence of furrows or deviations from a smooth curve
%by comparing the convex hull area to the object area. Clumped cells can be
%thought of a overlapping ellipses that would produce object furrows. It is
%a measure that is 0<x<1 and has a similar distribution to the form-factor.
[sfarray,sfthresh] = solidityfilter(nucleiSolidity_array);
myplothistformfactor(nucleiSolidity_array,'nucleiSolidity',sfthresh)
%Compactness gives a second way to measure furrows or spikes. Compactness
%should be greater than 1. Taking the log of this data will move the
%smallest value of the disribution to 0 and make it easier to compare
%outliers, which might have extreme values since this statistic is a ratio.
[cfarray,cfthresh] = compactnessfilter(nucleiCompactness_array);
myplothistformfactor(nucleiCompactness_array,'nucleiCompactness',cfthresh)
nucleiLogic = afarray & cfarray & fffarray & sfarray;
myData = nucleiMeanIntensityDAPI_array.*nucleiArea_array;
myData = myData(nucleiLogic);
myplothistDAPI(myData,'TotalNucleiDAPI')
myData = nucleiMeanIntensityYFP_array.*nucleiArea_array;
myData = myData(nucleiLogic);
myplothist(myData,'TotalNucleiYFP')
myData = nucleiMeanIntensityCy5_array.*nucleiArea_array;
myData = myData(nucleiLogic);
myplothist(myData,'TotalNucleiCy5')

my_plotregression(nucleiMeanIntensityYFP_array,nucleiMeanIntensityCy5_array)
my_plotcoefofvar(nucleiMeanIntensityYFP_array,nucleiMeanIntensityCy5_array)
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

function [] = my_plotregression(nucleiMeanIntensityYFP_array,nucleiMeanIntensityCy5_array)
h=figure ( 'visible', 'off', 'position', [10, 10, 672, 512] );
ax = plotregression(log(nucleiMeanIntensityYFP_array'),log(nucleiMeanIntensityCy5_array'));
xlabel(ax,'p53-YFP (Target)')
saveas (ax, 'linearregression', 'pdf' );
close(h);
end

function [] = my_plotcoefofvar(nucleiMeanIntensityYFP_array,nucleiMeanIntensityCy5_array)
h=figure ( 'visible', 'off', 'position', [10, 10, 672, 512] );
ax=axes('parent',h);
datain = log(nucleiMeanIntensityCy5_array);
[~,yyy] = hist(datain,20);
p53axis = zeros(1,19);
p53cov = zeros(1,19);
for i = 2:20
    temp = log(nucleiMeanIntensityYFP_array(((datain<yyy(i)) & (datain<yyy(i-1)))));
    if length(temp)<20
        p53axis(i-1) = NaN;
        p53cov(i-1) = NaN;
    else
        p53std = std(temp);
        p53m = mean(temp);
        p53axis(i-1) = (yyy(i)-yyy(i-1))/2+yyy(i-1);
        p53cov(i-1) = abs(p53std/p53m);
    end
end
plot(p53axis,p53cov,'.','MarkerSize',20)
xlabel(ax,'log(p53-YFP) A.U.')
ylabel(ax,'p53-YFP coefficient of variation')
saveas (ax, 'cov', 'pdf' );
close(h);
end

function [] = myplothist(in,filename)
h=figure ( 'visible', 'off', 'position', [10, 10, 672, 512] );
ax=axes('parent',h);
[y,x2] = hist(in,100);
h2 = bar(ax,x2,y/(sum(y)*(x2(2)-x2(1))),'hist');
set(h2,'FaceColor',[0 0.5 1],'EdgeColor',[0 0 1]);
title(filename);
hold on
y2 = ksdensity(in(in>0),x2,'support','positive');
plot(ax,x2,y2,'Color','black','LineWidth',1.5)
hold off
saveas (ax, filename, 'pdf' );
close(h);
end

function [] = myplothistDAPI(in,filename)
%Trim the DAPI data to represent FACs data
[y,xi] = ksdensity(in,'support','positive');
threshRight = xi(triangleThreshCore(y));
[y,xi] = ksdensity(-in);
threshLeft = -xi(triangleThreshCore(y));
inlogic = in>threshLeft & in<threshRight;
in = in(inlogic);
%Plot the data
h=figure ( 'visible', 'off', 'position', [10, 10, 672, 512] );
ax=axes('parent',h);
[y,x2] = hist(in,100);
h2 = bar(ax,x2,y/(sum(y)*(x2(2)-x2(1))),'hist');
set(h2,'FaceColor',[0 0.5 1],'EdgeColor',[0 0 1]);
title(filename);
hold on
y2 = ksdensity(in(in>0),x2,'support','positive','kernel','triangle');
plot(ax,x2,y2,'Color','black','LineWidth',1.5)
hold off
saveas (ax, filename, 'pdf' );
close(h);
end

function [] = myplothistarea(x,filename,vline,mu,sigma)
h=figure ( 'visible', 'off', 'position', [10, 10, 672, 512] );
ax=axes('parent',h);
[y,x2] = hist(x,100);
h2 = bar(ax,x2,y/(sum(y)*(x2(2)-x2(1))),'hist');
set(h2,'FaceColor',[0 0.5 1],'EdgeColor',[0 0 1]);
title(filename);
hold on
y2 = lognpdf(x2,mu,sigma);
plot(ax,x2,y2,'Color','black','LineWidth',1.5)
yL = get(ax,'YLim');
line([vline vline],yL,'Color','r','LineWidth',1.5);
hold off
saveas (ax, filename, 'pdf' );
close(h);
end

function [] = myplothistformfactor(in,filename,vline)
h=figure ( 'visible', 'off', 'position', [10, 10, 672, 512] );
ax=axes('parent',h);
[y,x2] = hist(in,100);
h2 = bar(ax,x2,y/(sum(y)*(x2(2)-x2(1))),'hist');
set(h2,'FaceColor',[0 0.5 1],'EdgeColor',[0 0 1]);
title(filename);
hold on
y2 = ksdensity(in(in>0),x2,'support','positive');
plot(ax,x2,y2,'Color','black','LineWidth',1.5)
yL = get(ax,'YLim');
line([vline(1) vline(1)],yL,'Color','r','LineWidth',1.5);
line([vline(2) vline(2)],yL,'Color','r','LineWidth',1.5);
hold off
saveas (ax, filename, 'pdf' );
close(h);
end

function [afarray,afthresh,meanlogx,sigmalogx] = areafilter(in)
L = length(in);
%assume area is distributed as a log-normal distribution
logx = log(in);
%assume there are outliers. We want to ignore these outliers while fitting
%the normal distribution, so we discard the bottom and top 10%.
logx = sort(logx);
logx = logx(round(L*.1):round(L*.9));
%From this trimmed set we have to bootstrap the tails of the distribution.
bootset = randn(L,1);
bootset = sort(bootset);
%We must estimate the variance of the distribution from the trimmed set.
%This variance must be scaled by the variance of a zero mean, unit variance
%normal distribution trimmed the same way. The z-score for a 10% tail is
%1.2815.
%The variance for this trimmed distribution can be found numerically using
%MATLAB. x = randn(1e6,1);x = sort(x);x = x(round(length(x)*.1):round(length(x)*.9));var(x)
%he variance equals 0.4377
varlogx = var(logx)/.4377;
sigmalogx = sqrt(varlogx);
meanlogx = mean(logx);
%The bootset is scaled by the mean and variance of the sample data
bootset = bootset*sigmalogx+meanlogx;
%The bootstraped tails are added to the trimmed distribution.
bootmin = bootset(bootset<min(logx));
bootmax = bootset(bootset>max(logx));
newdist = [bootmin',logx,bootmax'];
%The new distribution is scaled to zero mean and unit variance
newdistmean = mean(newdist);
newdistvar = var(newdist);
newdist = (newdist-newdistmean)/sqrt(newdistvar);
%The goodness of fit to a normal distribution is tested by the
%kolmogorov-smirnoff test.
goodfitbool = kstest(newdist);
if goodfitbool
    disp('The distribution of nuclei areas is not log-normal. Check to see if it is bi-modal.')
end
%Use 3 sigma below the mean as the cutoff for nucleus size
afthresh = exp(meanlogx - 3*sqrt(varlogx));
afarray = in>afthresh;
end

function [fffarray,fffthresh] = formfactorfilter(in)
x = -log(in);
x = x(x>0);
[y,xi] = ksdensity(x,'support','positive');
fffthresh(2) = exp(0);
fffthresh(1) = exp(-xi(triangleThreshCore(y)));
fffarray = in>fffthresh(1) & in<fffthresh(2);
end

function [cfarray,cfthresh] = compactnessfilter(in)
x = in(in>1);
x = log(x);
[y,xi] = ksdensity(x,'support','positive');
cfthresh(1) = exp(0);
cfthresh(2) = exp(xi(triangleThreshCore(y)));
cfarray = in>cfthresh(1) & in<cfthresh(2);
end

function [sfarray,sfthresh] = solidityfilter(in)
x = -log(in);
x = x(x>0);
[y,xi] = ksdensity(x,'support','positive');
thresh1 = xi(triangleThreshCore(y));
x2 = x(x>thresh1);
[y,xi] = ksdensity(x2,'support','positive');
sfthresh(2) = exp(0);
sfthresh(1) = exp(-xi(triangleThreshCore(y)));
sfarray = in>sfthresh(1) & in<sfthresh(2);
end

function [ind2] = triangleThreshCore(n)
%Find the highest peak the histogram
[c,ind]=max(n);
%Assume the long tail is to the right of the peak and envision a line from
%the top of this peak to the end of the histogram.
%The slope of this line, the hypotenuse, is calculated.
x1=0;
y1=c;
x2=length(n)-ind;
y2=n(end);
m=(y2-y1)/(x2-x1); %The slope of the line

%----- Find the greatest distance -----
%We are looking for the greatest distance betweent the histrogram and line
%of the triangle via perpendicular lines
%The slope of all lines perpendicular to the histogram hypotenuse is the
%negative reciprocal
p=-1/m; %The slope is now the negative reciprocal
%We now have two slopes and two points for two lines. We now need to solve
%this two-equation system to find their intersection, which can then be
%used to calculate the distances
iarray=(0:(length(n)-ind));
L=zeros(size(n));
for i=iarray
    intersect=(1/(m-p))*[-p,m;-1,1]*[c;n(i+ind)-p*i];
    %intersect(1)= y coordinate, intersect(2)= x coordinate
    L(i+ind)=sqrt((intersect(2)-i)^2+(intersect(1)-n(i+ind))^2);
end
[~,ind2]=max(L);
end