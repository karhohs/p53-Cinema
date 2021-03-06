function [] = smfishStackPreprocessing(stackpath,varargin)
%noise reduction by matched filtering with a guassian kernel
%create gaussian filter that approximates 3D PSF of the microscope.
%Distance units are in micrometers.
%----- Set Parameters -----
parameters.objective = 60; %as in 60x
parameters.NA = 1.4; %typical of 60x oil immersion objections
parameters.rindex = 1.51; %typical refractive index of oil
parameters.camerapixellength = 6.45; %Both cameras in the Lahav have pixel dimensions of 6.45 x 6.45 um.
parameters.zstepsize = 0.3; %User defined with z-stack is obtained
parameters.wavelength = .67; %Cy5 probe wavelength approximately 670 nanometers

%----- Parse varargin -----
%'flatfieldpath' = 'path\2\files'
%'subtractbackground' = true
p = inputParser;
p.addRequired('stackpath', @(x)ischar(x));
p.addParamValue('flatfieldpath', '', @(x)ischar(x));
p.addParamValue('subtractbackground', true, @(x)islogical(x));
p.addParamValue('fluorchan','Cy5',@(x)ischar(x));
p.parse(stackpath, varargin{:});
%----- flatfield correction for each z-slice -----
if ~isempty(p.Results.flatfieldpath)
    flatfieldcorrection(stackpath,p.Results.flatfieldpath)
    tempfoldername=regexp(stackpath,'(?<=\\)[\w ]*','match'); %Prepare to create a new folder to place background subtracted stacks
    tempfoldername=[tempfoldername{end},'_ff'];
    stackpath=[stackpath,'\..\',tempfoldername];
end
%----- read the contents of the input folder -----
disp(['working in ', stackpath]); %Sanity Check
dirCon_stack = dir(stackpath);
%----- Create a new folder to hold corrected images -----
tempfoldername=regexp(stackpath,'(?<=\\)[\w ]*','match'); %Prepare to create a new folder to place background subtracted stacks
tempfoldername=[tempfoldername{end},'_smFISH'];
smfishstackpath=[stackpath,'\..\',tempfoldername];
mkdir(smfishstackpath);
cd(smfishstackpath)
smfishstackpath=pwd;
%process only the stacks that contain single molecule FISH
stacknames=importStackNames(dirCon_stack,p.Results.fluorchan);
if isempty(stacknames)
    disp('Did not find any z-stacks of the specified fluorescent channel.')
    return
end
stacknames2=stacknames; %I apologize for the really confusing file naming system.
%This is part of a bug. This file expects the input filename to be of a
%certain format. The following for-loop partially ensures the file names
%are in that format.
for i=1:length(stacknames)
    Name_temp = regexprep(stacknames(i),'\s',''); %Remove all not(alphabetic, numeric, or underscore) characters
    Name_temp = regexprep(Name_temp,'tocamera','','ignorecase'); %remove 'tocamera' if present b/c it is not informative
    Name_temp = regexprep(Name_temp,'camera','','ignorecase'); %remove 'camera' if present b/c it is not informative
    stacknames2(i) = Name_temp;
end

for bigInd = 1:length(stacknames)
    %----- Load the image file -----
    parameters.stacknametest = [stackpath '\' stacknames{bigInd}];
    [IM,sizeOfImage,hLoG,tempI1,tempI2,hMeanxy,hMeanz,IMMeanIntensity,hGaus,xy,z,pixelRatio] = variableInitialization(parameters);
    IM = loadZstack([stackpath '\' stacknames{bigInd}],IM,sizeOfImage);
    dataName = regexprep(stacknames2{bigInd},'(?<=_t)(\w*)(?=\.)','$1_data');
    dataName = regexp(dataName,'.*(?=\.)','match','once');
    save(dataName,'sizeOfImage');
    %----- background subtraction for each z-slice -----
    %I've heard that 3D deconvolution using an iterative blind-maximum-likehood
    %algorithm is very effective at removing out of focus light from each
    %z-slice. I played around with the AutoQuant software package at the NIC.
    %AutoQuant has this deconvolution algorithm and can batch process a whole
    %folder of files. The results were quite nice. However, the deconvolution
    %process is time consuming: it takes about 20 minutes for a 1344x1024x50
    %TIFF file and requires scheduling time at the NIC and using a computer at
    %the NIC for the processing. Perhaps later I can incorporate deconv. into
    %this MATLAB pipeline. Until then, as an alternative, I have found that
    %using the tried and true 'imopen' for background subtraction does a pretty
    %good job at removing out of focus light as long as the structuring element
    %is of an appropriate size (i.e. slightly bigger than a diffraction limited
    %spot).
    if p.Results.subtractbackground
        IM = JaredsBackground(IM);
    end
    %----- enhance diffraction limited spots using the LoG filter -----
    IM = imfilter(IM,hLoG,'symmetric'); %Totally works. sweet!
    %----- calculate the test statistics that will identify legit spots -----
    %%%Test Statistic 1: The 3D curvature. Gives a sense about how much a spot
    %resembles a point source of light. It gives a sense of the spots geometry
    %as opposed to the brightness of the spot.
    curvature = mySobelHessianCurvature(IM,tempI1,tempI2,pixelRatio);
    %A really large negative value indicates geometry like a point source. The
    %numbers produced are often extremely large and it may be a good idea to
    %normalize.
    %I turned a negative into a positive; I want you all to know that.
    curvature = -curvature;
    curvature(curvature<0) = 0;
    %%%Test Statistic 2: The mean brightness of the area. Indeed, we expect the
    %mRNA FISH signal to be brighter than the background. Taking the mean
    %reduces the weight of random peaks due to noise, since noise in these
    %images is of the zero-mean variety.
    %find local maxima
    fociCandidates = imregionalmax(IM,26);
    %Find the mean of a local volume that will capture an entire point source.
    for i=1:sizeOfImage(2)
        tempI1(:,i,:) = imfilter(reshape(IM(:,i,:),[sizeOfImage(1), sizeOfImage(3)]),hMeanz,'symmetric'); %z
    end
    for i=1:sizeOfImage(3)
        tempI2(:,:,i) = imfilter(tempI1(:,:,i),hMeanxy,'symmetric'); %y
    end
    for i=1:sizeOfImage(3)
        IMMeanIntensity(:,:,i) = imfilter(tempI2(:,:,i),hMeanxy','symmetric'); %x
    end
    %The final test statistic is the product of test statistic 1 and 2
    spotStat = IMMeanIntensity.*curvature;
    %----- Find a threshold that separates signal from noise -----
    spotStat = spotStat.*fociCandidates;
    index = find(spotStat);
    %clear fociCandidates
    if iscolumn(index)
        index = index';
    end
    spotStat2 = spotStat(index);
    IMMeanIntensity2 = IMMeanIntensity(index);
    curvature2 = curvature(index);
    %The test statistic tends to vary over several orders of magnitude
    %therefore it is easier to compare these values in a log scale.
    spotStat2 = log(spotStat2);
    curvature2 = log(curvature2);
    if isrow(spotStat2)
        spotStat2 = spotStat2';
    end
    if isrow(curvature2)
        curvature2 = curvature2';
    end
    if isrow(IMMeanIntensity2)
        IMMeanIntensity2 = IMMeanIntensity2';
    end
    if isrow(index)
        index = index';
    end
    spotStat3 = [spotStat2,index,curvature2,IMMeanIntensity2];
    spotStat3 = sortrows(spotStat3,1);
    save(dataName,'spotStat3','-append');
    %find a good threshold
    [threshold,n,xout,n2,xout2] = triminthresh(spotStat3(:,1)); %#ok<NASGU,ASGLU>
    save(dataName,'threshold','n','xout','n2','xout2','-append');
    ind = find(spotStat3(:,1)>threshold,1,'first');
    foci = zeros(sizeOfImage);
    foci(spotStat3(ind:end,2)) = 1;
    fociarray = spotStat3(ind:end,2);
    save(dataName,'fociarray','-append');
    %----- Create the final image with bonafide spots and other aesthetic images -----
    %sum projection of foci
    sumProj = sum(foci,3);
    Name = regexprep(stacknames2{bigInd},'(\w*)(?=\.)','$1_sumProj');
    imwrite(uint8(sumProj),[smfishstackpath,'\',Name],'tif','WriteMode','append','Compression','none');
    %max project the stamp
    stampProj3D = padarray(zeros(sizeOfImage),[xy xy z],'symmetric');
    %gaussian stamp (approx. the PSF) on the sum projection of foci
    foci2 = padarray(foci,[xy xy z]);
    index = find(foci2);
    if iscolumn(index)
        index = index';
    end
    s = size(foci2);
    for i=index
        [i2,j2,k2] = ind2sub(s,i);
        stampProj3D(i2-xy:i2+xy,j2-xy:j2+xy,k2-z:k2+z) = stampProj3D(i2-xy:i2+xy,j2-xy:j2+xy,k2-z:k2+z) + hGaus;
    end
    stampProj3D2 = stampProj3D(xy+1:end-xy,xy+1:end-xy,z+1:end-z);
    %Save the 3D image
    %     Name = regexprep(stacknames2{bigInd},'(?<=_t)(\w*)(?=\.)','$1_stampProj3D');
    %     for i = 1:sizeOfImage(3)
    %         imwrite(uint8(stampProj(:,:,i)),[smfishstackpath,'\',Name],'tif','WriteMode','append','Compression','none');
    %     end
    %Project the 3D image into 2D
    stampProj = sum(stampProj3D2,3);
    Name = regexprep(stacknames2{bigInd},'(\w*)(?=\.)','$1_stampProj');
    imwrite(uint8(stampProj),[smfishstackpath,'\',Name],'tif','WriteMode','append','Compression','none');
    
    %Create Max Projection of the input image
    maxProj = max(IM,[],3);
    maxProj = uint16(maxProj);
    Name = regexprep(stacknames2{bigInd},'(\w*)(?=\.)','$1_maxProj');
    imwrite(uint16(maxProj),[smfishstackpath,'\',Name],'tif','WriteMode','append','Compression','none');
    %Create Merged Color image
    maxProj = bitshift(maxProj, -4);
    maxProj = uint8(maxProj);
    [s1 s2] = size(maxProj);
    maxProj2 = zeros(s1,s2,3);
    maxProj2(:,:,1) = maxProj;
    maxProj2(:,:,2) = maxProj;
    maxProj2(:,:,3) = maxProj;
    for i = 1:length(fociarray)
        [y2,x2,~] = ind2sub(sizeOfImage,fociarray(i));
        maxProj2(y2,x2,:) = [255 0 0];
    end
    Name = regexprep(stacknames2{bigInd},'(\w*)(?=\.)','$1_ColorMerge');
    imwrite(uint8(maxProj2),[smfishstackpath,'\',Name],'tif','WriteMode','append','Compression','none');
    %3D scatter plot
    %[y2,x2,z2] = ind2sub(s,fociarray);
    %scatter3(x2,y2,z2)
    Name = regexprep(stacknames2{bigInd},'(\w*)(?=\.)','$1_maxProj');
    smfishPlot([dataName '.mat'],smfishstackpath,Name,stacknames2{bigInd});
end
signalCompletionWithEmail();
signalCompletionWithSound();
end

function [tempI1] = mySobelHessianCurvature(I,tempI1,tempI2,pixelRatio)
%This function uses the Sobel filter to approximate the derivatives in a
%gradient. Since the sobel filter is seperable the it can also be
%conveniently extended to find the Hessian.
%The image I is 3D and has coordinates (y,x,z).
h1 = [0.25 0.5 0.25];
h2 = [-0.5 0 0.5];
h2z = h2/pixelRatio;
[sy,sx,sz] = size(I);
Fx = zeros(size(I));
Fy = zeros(size(I));
Fz = zeros(size(I));
Fxx = zeros(size(I));
Fxy = zeros(size(I));
Fxz = zeros(size(I));
Fyy = zeros(size(I));
Fyz = zeros(size(I));
Fzz = zeros(size(I));
%The Sobel separted filters to find the Fx
for i=1:sx
    tempI1(:,i,:) = imfilter(reshape(I(:,i,:),[sy sz]),h1,'symmetric'); %z
end
for i=1:sz
    tempI2(:,:,i) = imfilter(tempI1(:,:,i),h1,'symmetric'); %y
end
for i=1:sz
    Fx(:,:,i) = imfilter(tempI2(:,:,i),h2','symmetric'); %x
end
%The Sobel separted filters to find the Fy
for i=1:sx
    tempI1(:,i,:) = imfilter(reshape(I(:,i,:),[sy sz]),h1,'symmetric'); %z
end
for i=1:sz
    tempI2(:,:,i) = imfilter(tempI1(:,:,i),h2,'symmetric'); %y
end
for i=1:sz
    Fy(:,:,i) = imfilter(tempI2(:,:,i),h1','symmetric'); %x
end
%The Sobel separted filters to find the Fz
for i=1:sx
    tempI1(:,i,:) = imfilter(reshape(I(:,i,:),[sy sz]),h2z,'symmetric'); %z
end
for i=1:sz
    tempI2(:,:,i) = imfilter(tempI1(:,:,i),h1,'symmetric'); %y
end
for i=1:sz
    Fz(:,:,i) = imfilter(tempI2(:,:,i),h1','symmetric'); %x
end
%The Sobel separted filters to find the Fxx
for i=1:sx
    tempI1(:,i,:) = imfilter(reshape(Fx(:,i,:),[sy sz]),h1,'symmetric'); %z
end
for i=1:sz
    tempI2(:,:,i) = imfilter(tempI1(:,:,i),h1,'symmetric'); %y
end
for i=1:sz
    Fxx(:,:,i) = imfilter(tempI2(:,:,i),h2','symmetric'); %x
end
%The Sobel separted filters to find the Fxy
for i=1:sx
    tempI1(:,i,:) = imfilter(reshape(Fx(:,i,:),[sy sz]),h1,'symmetric'); %z
end
for i=1:sz
    tempI2(:,:,i) = imfilter(tempI1(:,:,i),h2,'symmetric'); %y
end
for i=1:sz
    Fxy(:,:,i) = imfilter(tempI2(:,:,i),h1','symmetric'); %x
end
%The Sobel separted filters to find the Fxz
for i=1:sx
    tempI1(:,i,:) = imfilter(reshape(Fx(:,i,:),[sy sz]),h2z,'symmetric'); %z
end
for i=1:sz
    tempI2(:,:,i) = imfilter(tempI1(:,:,i),h1,'symmetric'); %y
end
for i=1:sz
    Fxz(:,:,i) = imfilter(tempI2(:,:,i),h1','symmetric'); %x
end
%The Sobel separted filters to find the Fyy
for i=1:sx
    tempI1(:,i,:) = imfilter(reshape(Fy(:,i,:),[sy sz]),h1,'symmetric'); %z
end
for i=1:sz
    tempI2(:,:,i) = imfilter(tempI1(:,:,i),h2,'symmetric'); %y
end
for i=1:sz
    Fyy(:,:,i) = imfilter(tempI2(:,:,i),h1','symmetric'); %x
end
%The Sobel separted filters to find the Fyz
for i=1:sx
    tempI1(:,i,:) = imfilter(reshape(Fy(:,i,:),[sy sz]),h2z,'symmetric'); %z
end
for i=1:sz
    tempI2(:,:,i) = imfilter(tempI1(:,:,i),h1,'symmetric'); %y
end
for i=1:sz
    Fyz(:,:,i) = imfilter(tempI2(:,:,i),h1','symmetric'); %x
end
%The Sobel separted filters to find the Fzz
for i=1:sx
    tempI1(:,i,:) = imfilter(reshape(Fz(:,i,:),[sy sz]),h2z,'symmetric'); %z
end
for i=1:sz
    tempI2(:,:,i) = imfilter(tempI1(:,:,i),h1,'symmetric'); %y
end
for i=1:sz
    Fzz(:,:,i) = imfilter(tempI2(:,:,i),h1','symmetric'); %x
end
%Find the curvature matrix
myHessian = zeros(3);
for i = 1:numel(I)
    myHessian(1,1) = Fxx(i);
    myHessian(1,2) = Fxy(i);
    myHessian(1,3) = Fxz(i);
    myHessian(2,1) = Fxy(i);
    myHessian(2,2) = Fyy(i);
    myHessian(2,3) = Fyz(i);
    myHessian(3,1) = Fxz(i);
    myHessian(3,2) = Fyz(i);
    myHessian(3,3) = Fzz(i);
    tempI1(i) = det(myHessian);
end
end

function [threshold,n,xout,n2,xout2] = triminthresh(A)
%Calculate a few rank statistics (assumes A is already sorted)
la = length(A);
q1a = A(round(0.25*la)); %first quartile
q2a = A(round(0.50*la));
q3a = A(round(0.75*la)); %third quartile
myIQRa = q3a-q1a;
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
for i = 2:length(n2der);
    if (n2der(i-1)<0 && n2der(i)>=0)
        ind = i-1;
        break
    elseif (abs(n2der(i-1))<=1) && (n2(i-1) == 0 || n2(i-1) == 1 || n2(i-1) == 2)
        ind = i-1;
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

function [S] = JaredsBackground(S)
resizeMultiplier = 1/2; % Downsampling scale factor makes image processing go faster and smooths image
seSize2 = 40; % I find the value of 25 works well with 60x, binning 1, mRNA FISH images
se2 = strel('disk', seSize2*resizeMultiplier);  %Structing elements are necessary for using MATLABS image processing functions
origSize  = size(S);
for k=1:origSize(3)
    % Rescale image and compute background using closing/opening.
    I    = imresize(S(:,:,k), resizeMultiplier);
    pad   = ceil(seSize2*resizeMultiplier);
    % Pad image with a reflection so that borders don't introduce artifacts
    I    = padarray(I, [pad,pad], 'symmetric', 'both');
    % Perform opening/closing to get background
    I     = imopen(I, se2);     % ignore high-intensity features typical of mRNA spots
    % Remove padding and resize
    I     = floor(imresize(I(pad+1:end-pad, pad+1:end-pad), origSize(1:2)));
    S(:,:,k) = S(:,:,k) - I;
end
S(S<0)=0;
end

function [IM] = loadZstack(path,IM,s)
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

function [Temp] = importStackNames(dirCon_stack,fc)
expr=['.*(?<!thumb.*)_w\d+' fc '.*'];
Temp=cell([1,length(dirCon_stack)]); %Initialize cell array
% ----- Identify the legitimate stacks -----
i=1;
for j=1:length(dirCon_stack)
    Temp2=regexp(dirCon_stack(j).name,expr,'match','once','ignorecase');
    if Temp2
        Temp{i}=Temp2;
        i=i+1;
    end
end
% ----- Remove empty cells -----
Temp(i:end)=[];
% for j=length(Temp):-1:1
%     if isempty(Temp{j})
%         Temp(j)=[];
%     end
% end
end

function [IM,sizeOfImage,hLoG,tempI1,tempI2,hMeanxy,hMeanz,IMMeanIntensity,hGaus,xy,z,pixelRatio] = variableInitialization(p)
%This function was made in a effort to speed things up. I'm not sure it did
%that. It may have just made the code more difficult to read.
info = imfinfo(p.stacknametest,'tif');
sizeOfImage = [info(1).Height, info(1).Width, length(info)];

IM = zeros(sizeOfImage);
tempI1 = zeros(sizeOfImage);
tempI2 = zeros(sizeOfImage);
IMMeanIntensity = zeros(sizeOfImage);
sigmaXYos = 0.21*p.wavelength/p.NA; %lateral st. dev of the gaussian filter in object space
sigmaZos = 0.66*p.wavelength*p.rindex/(p.NA^2) ;%axial st. dev of the gaussian filter in object space
Pxy = p.camerapixellength/p.objective; %lateral pixel size
sigmaXY = sigmaXYos/Pxy; %lateral st. dev of gaussian filter in image space
sigmaZ = sigmaZos/p.zstepsize; %axial st. dev of gaussian filter in image space
xy = round(3*sigmaXY);
z = round(3*sigmaZ);
xyMLV = round(4*sigmaXY);
zMLV = round(4*sigmaZ);
K = 1/((2*pi)^(3/2)*sqrt(sigmaXY^2*sigmaZ)); % log3d coefficient
log3d = @(x,y,z) K*exp(-0.5*(x^2/sigmaXY+y^2/sigmaXY+z^2/sigmaZ))*((x^2-4*sigmaXY)/(4*sigmaXY^2)+(y^2-4*sigmaXY)/(4*sigmaXY^2)+(z^2-4*sigmaZ)/(4*sigmaZ^2));
hLoG = zeros(2*xy+1,2*xy+1,2*z+1);
for i=1:2*xy+1
    for j=1:2*xy+1
        for k=1:2*z+1
            hLoG(i,j,k) = log3d(i-1-xy,j-1-xy,k-1-z); %the 3D filter
        end
    end
end
%tune the filter so that it does not amplify the signal.
temp1 = ones(2*xy+1,2*xy+1,2*z+1);
hLoG=-hLoG; %otherwise the center weight, the largest weight, is negative.
temp2 = sum(sum(sum(temp1.*hLoG)));
hLoG=hLoG/temp2;
K2=1/(xyMLV*xyMLV*zMLV);
hMeanxy=ones(1,xyMLV);
hMeanz=K2*ones(1,zMLV);
%3D Gaussian Filter\Stamp\PSF approximation
mu = [0,0,0]; %zero mean gaussian
SIGMA = [sigmaXY,0,0;0,sigmaXY,0;0,0,sigmaZ];
hGaus = zeros(2*xy+1,2*xy+1,2*z+1);
for i=1:2*xy+1
    for j=1:2*xy+1
        for k=1:2*z+1
            hGaus(i,j,k) = mvnpdf([i-1-xy,j-1-xy,k-1-z],mu,SIGMA); %the 3D filter
        end
    end
end
hGaus=hGaus*(2^5)/(max(max(max(hGaus))));
pixelRatio = sigmaZ/sigmaXY;
end

function [] = signalCompletionWithSound()
global playerkwk
if ~isempty(playerkwk)
    play(playerkwk);
else
    disp('Victory over Data!')
end
end

function [] = signalCompletionWithEmail()
% Send the email
sendmail('gandalfisarockstar@gmail.com', 'Mail from MATLAB', ...
    'Hi Kyle! Your MATLAB run is complete!')
% sendmail('gandalfisarockstar@gmail.com', 'Mail from MATLAB', ...
%     'Hi Kyle! Your MATLAB run is complete!', ...
%     {'sub_folder/signals.m', 'system.mdl'})
end