function [] = lowCopyNumberThresholdFinder(smfishpath,varargin)
% [] = lowCopyNumberThresholdFinder(smfishpath,'manualthresh',false)
% Input:
% smfishpath: a char. The path to the folder that contains the output from
% smfishStackPreprocessing().
% manualthresh: a logical variable that when set to true will launch a
% user-interface that guides the selection of a threshold.
%
% Output:
% There is no direct argument output. Rather, a new folder will be created
% that recomputes the output of the smfishStackPreprocessing() function
% with a universal threshold that was calculated using all of the data in
% from the previously created output folder.
%
% Description:
% The idea is that the automatic detection of foci in the case where the
% number of foci per cell is very low, <5, can be improved by considering a
% whole collection of images instead of just one (that may only have a few
% cells to begin with). I have observed that many of the
% foci identified in low copy number situations are questionable, because
% they are either very dim or located extra-cellularly.
%
% There will also be an option to choose the threshold manually.
%
% Other Notes:
%

p = inputParser;
p.addRequired('smfishpath', @(x)ischar(x));
p.addParamValue('manualthresh', false, @(x)islogical(x));
p.parse(smfishpath, varargin{:});
if p.Results.manualthresh
    
end

%Read the directory where all the .mat files of fish data are kept.
%----- read the contents of the input folder -----
disp(['working in ', smfishpath]); %Sanity Check
dirCon_stack = dir(smfishpath);

%Create a new folder using the source folder name.
%----- Create a new folder to hold corrected images -----
tempfoldername=regexp(smfishpath,'(?<=\\)[\w ]*','match'); %Prepare to create a new folder to place background subtracted stacks
tempfoldername=[tempfoldername{end},'_lowCopyNum'];
smfishstackpath=[smfishpath,'\..\',tempfoldername];
mkdir(smfishstackpath);
cd(smfishstackpath)
% The .mat files can be identified by having "_data.mat" at the end of the filename.
expr='.*_data.mat';
datanames=cell([1,length(dirCon_stack)]); %Initialize cell array
% ----- Identify the legitimate stacks -----
i=1;
for j=1:length(dirCon_stack)
    Temp2=regexp(dirCon_stack(j).name,expr,'match','once','ignorecase');
    if Temp2
        datanames{i}=Temp2;
        i=i+1;
    end
end
% ----- Remove empty cells -----
datanames(i:end)=[];
% Read each .mat file and store each "spotStat3" variable in a cell.
% spotStat3 has four columns organized as follows:
%[log(intensity*curvature), index, curvature, intensity];
%We are interested in creating an array with all of the
%log(intensity*curvature) data. Create that array
spotStatCell = cell(size(datanames));
numfoci_before = 0;
for i=1:length(spotStatCell)
    spotStatTemp = load(fullfile(smfishpath,datanames{i}), 'spotStat3');
    spotStatCell{i} = spotStatTemp.spotStat3(:,1);
    spotStatTemp2 = load(fullfile(smfishpath,datanames{i}), 'fociarray');
    numfoci_before = numfoci_before + length(spotStatTemp2.fociarray);
end
spotStatArray = cell2mat(spotStatCell');
%Use the triminthresh() method to calculate a universal threshold.
[threshold,~,~,~,~] = triminthresh(spotStatArray);
%Count foci and create images using the new threshold
numfoci_after = 0;
for i = 1:length(datanames)
    spotStatTemp = load(fullfile(smfishpath,datanames{i}), 'spotStat3');
    dataName = datanames{i};
    spotStat3 = spotStatTemp.spotStat3; %#ok<NASGU>
    save(dataName,'spotStat3');
    %identify foci using the universal threshold
    sizeOfImage = load(fullfile(smfishpath,datanames{i}), 'sizeOfImage');
    sizeOfImage = sizeOfImage.sizeOfImage;
    save(dataName,'sizeOfImage','-append');
    save(dataName,'threshold','-append');
    ind = find(spotStatTemp.spotStat3(:,1)>threshold,1,'first');
    foci = zeros(sizeOfImage);
    foci(spotStatTemp.spotStat3(ind:end,2)) = 1;
    fociarray = spotStatTemp.spotStat3(ind:end,2);
    numfoci_after = numfoci_after + length(fociarray);
    save(dataName,'fociarray','-append');
    %sum projection of foci
    sumProj = sum(foci,3);
    Name = regexprep(datanames{i},'(\w*)_data.*','$1_sumProj.tif');
    imwrite(uint8(sumProj),[smfishstackpath,'\',Name],'tif','WriteMode','append','Compression','none');
    pic = rings(fociarray,sizeOfImage);
    Name = regexprep(datanames{i},'(\w*)_data.*','$1_circles.tif');
    imwrite(uint8(pic),[smfishstackpath,'\',Name],'tif','WriteMode','append','Compression','none');
end
str = sprintf('foci before = %d and foci after = %d \n The new threshold is %f', numfoci_before, numfoci_after, threshold);
disp(str)

%Generate a report
hist(spotStatArray,200)
title('global')
spotStatMean = zeros(size(spotStatCell));
for i=1:length(datanames)
    spotStatTemp = load(fullfile(smfishpath,datanames{i}), 'threshold');
    thresh = spotStatTemp.threshold;
    spotStatMean(i) = mean(spotStatCell{i}(spotStatCell{i}>thresh));
    figure
    hist(spotStatCell{i},200)
    [thresh2,~,~,~,~] = triminthresh(spotStatCell{i});
    str = sprintf('%d: thresh = %f, thresh2 = %f',i,thresh,thresh2);
    title(str)
end
figure
plot(spotStatMean,'o','MarkerEdgeColor','k','MarkerFaceColor',[.49 1 .63], 'MarkerSize',8)
end

function [threshold,n,xout,n2,xout2] = triminthresh(A)
A = sort(A);
la = length(A);
q1a = A(round(0.25*la)); %first quartile
q2a = A(round(0.50*la));
myIQRa = 2*(q2a-q1a);
myCutoffa = 3*myIQRa+q2a;
%Create the histogram
[n,xout]=hist(A,100);
%Use the triangle threshold for the initial guess
ind = triangleThreshCore(n);
threshold = xout(ind);
%Look for minimum change in the number of foci or when the change in foci
%is less than 1.
B = A(A>threshold);
[n2,xout2] = hist(B,100);
n2der = smooth(n2);
n2der = conv(n2der,[0.5 0 -0.5],'same'); %the central difference derivative to find the min
sumn2 = zeros(size(n2));
for i=1:length(n2)
    sumn2(i) = sum(n2(i:end));
end
p = n2der'./sumn2;
for i = 2:length(n2der);
    if (n2der(i-1)<0 && n2der(i)>=0) %tests for a minimum using the first derivative
        ind = i-1;
        break
    elseif (abs(n2der(i-1))<=1) && (n2(i-1) == 0 || n2(i-1) == 1 || n2(i-1) == 2) %tests for a flat region where the change in foci is less than 1
        ind = i-1;
        for j = 2:length(n2der)
            if (abs(p(j))<0.03 && ((p(j)-p(j-1))>0)) %tests for a minimal about of change: less than 3% of foci beyond hypothetical threshold
                ind2 = j-1;
                if ind2<ind
                    ind = ind2;
                end
                break
            end
        end
        break
        
    end
end
threshold = xout2(ind);
logicStepCounter = 1;
while logicStepCounter ~= 0
    switch logicStepCounter
        case 1
            if threshold<myCutoffa
                %Find the peak of the putative signal.
                [~,putativeSignalPeakInd] = max(n2);
                putativeSignalPeak = xout2(putativeSignalPeakInd);
                logicStepCounter = 2;
            else
                break
            end
        case 2
            if (putativeSignalPeak>myCutoffa) && (putativeSignalPeakInd>=2)
                %If the peak is greater than the cutoff than go forward with
                %the threshold search. Look for the min to the left of this peak.
                %It will be the final threshold.
                for i = putativeSignalPeakInd:-1:2;
                    if (n2der(i-1)<0 && n2der(i)>=0)
                        ind = i-1;
                        break
                    elseif (abs(n2der(i-1))<=1) && (n2(i-1) == 0 || n2(i-1) == 1 || n2(i-1) == 2)
                        ind = i-1;
                        break
                    end
                end
                threshold = xout2(ind);
                break
            else
                logicStepCounter = 3;
            end
        case 3
            %Repeat the triangle threshold method
            C = B(B>putativeSignalPeak);
            [n3,xout3] = hist(C,100);
            ind = triangleThreshCore(n3);
            threshold = xout3(ind);
            C = B(B>threshold);
            [n3,xout3] = hist(C,100);
            n3der = smooth(n3);
            n3der = conv(n3der,[0.5 0 -0.5],'same'); %the central difference derivative to find the min
            n3der = smooth(n3der);
            for i = 2:length(n3der);
                if (n3der(i-1)<0 && n3der(i)>=0)
                    ind = i-1;
                    break
                elseif (abs(n3der(i-1))<=1) && (n3(i-1) == 0 || n3(i-1) == 1 || n3(i-1) == 2)
                    ind = i-1;
                    break
                end
            end
            threshold = xout3(ind);
            break
        otherwise
            disp('If you see this message something went wrong during threshold calculation.');
            break
    end
end
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

function [ringPicture] = rings(fociarray,sizeOfImage)
linewidth = 1;
circleradiusbase = 5; %units in pixels
circlecircumference = 2*pi*circleradiusbase;
ringIM = zeros(sizeOfImage(1:2));
for k = 1:linewidth;
    if linewidth > 1
        circleradius = circleradiusbase + k - linewidth/2;
    else
        circleradius = circleradiusbase;
    end
    for i = 1:length(fociarray)
        [y2,x2,~] = ind2sub(sizeOfImage,fociarray(i));
        for j=1:round(2*circlecircumference)
            circlex = round(x2 + circleradius*cos(j*2*pi/round(2*circlecircumference)));
            circley = round(y2 + circleradius*sin(j*2*pi/round(2*circlecircumference)));
            if circlex < 1
                continue
            elseif circlex > sizeOfImage(2)
                continue
            end
            if circley < 1
                continue
            elseif circley > sizeOfImage(1)
                continue
            end
            ringIM(circley,circlex) = 1;
        end
    end
end
ringIM = logical(ringIM);
ringPicture = zeros(sizeOfImage(1),sizeOfImage(2),3);
ringPicture(:,:,1) = 255*ringIM;
ringPicture(:,:,2) = 255*ringIM;
end