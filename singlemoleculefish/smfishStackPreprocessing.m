function [] = smfishStackPreprocessing()
%noise reduction by matched filtering with a guassian kernel
%create gaussian filter that approximates 3D PSF of the microscope.
%Distance units are in micrometers.
%----- Set Parameters -----
objective = 60; %as in 60x
NA = 1.4; %typical of 60x oil immersion objections
rindex = 1.51; %typical refractive index of oil
camerapixellength = 6.45; %Both cameras in the Lahav have pixel dimensions of 6.45 x 6.45 um.
zstepsize = 0.25; %User defined with z-stack is obtained
wavelength = .67; %Cy5 probe wavelength approximately 670 nanometers
%----- Load the image file -----
path = 'M:\SI_201203B\deconvolution_test\orignal_data_ff_bkgd\20120312_pos1-1_w2Cy5_s1_t1_ff.TIF';
info = imfinfo(path,'tif');
t = Tiff(path,'r');
IM_temp = t.read;
[leny,lenx] = size(IM_temp);
lenz = length(info);
IM = zeros(leny,lenx,lenz);
if lenz > 1
    for k=1:length(info)-1
        IM(:,:,k) = double(t.read);
        t.nextDirectory;
    end
    %one last time without t.nextDirectory
    IM (:,:,lenz) = double(t.read);
end
t.close;
%----- flatfield correction for each z-slice -----
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

%----- enhance diffraction limited spots using the LoG filter -----
sigmaXYos = 0.21*wavelength/NA; %lateral st. dev of the gaussian filter in object space
sigmaZos = 0.66*wavelength*rindex/(NA^2) ;%axial st. dev of the gaussian filter in object space
Pxy = camerapixellength/objective; %lateral pixel size
sigmaXY = sigmaXYos/Pxy; %lateral st. dev of gaussian filter in image space
sigmaZ = sigmaZos/zstepsize; %axial st. dev of gaussian filter in image space
xy = round(3*sigmaXY);
z = round(3*sigmaZ);
K = 1/((2*pi)^(3/2)*sqrt(sigmaXY^2*sigmaZ)); % log3d coefficient
log3d = @(x,y,z) K*exp(-0.5*(x^2/sigmaXY+y^2/sigmaXY+z^2/sigmaZ))*((x^2-4*sigmaXY)/(4*sigmaXY^2)+(y^2-4*sigmaXY)/(4*sigmaXY^2)+(z^2-4*sigmaZ)/(4*sigmaZ^2));
h = zeros(2*xy+1,2*xy+1,2*z+1);
for i=1:2*xy+1
    for j=1:2*xy+1
        for k=1:2*z+1
            h(i,j,k) = log3d(i-1-xy,j-1-xy,k-1-z); %the 3D filter
        end
    end
end
%tune the filter so that it does not amplify the signal.
temp1 = ones(2*xy+1,2*xy+1,2*z+1);
h=-h; %otherwise the center weight, the largest weight, is negative.
temp2 = temp1.*h;
h=h/temp2;
clear temp1;
clear temp2;
IM = imfilter(IM,h,'symmetric'); %Totally works. sweet!
%----- calculate the test statistics that will identify legit spots -----
%%%Test Statistic 1: The 3D curvature. Gives a sense about how much a spot
%resembles a point source of light. It gives a sense of the spots geometry
%as opposed to the brightness of the spot.
curvature = mySobelHessianCurvature(IM);
%A really large negative value indicates geometry like a point source. The
%numbers produced are often extremely large and it may be a good idea to
%normalize.
%I turned a negative into a positive; I want you all to know that.
curvature = abs(curvature);
%%%Test Statistic 2: The mean brightness of the area. Indeed, we expect the
%mRNA FISH signal to be brighter than the background. Taking the mean
%reduces the weight of random peaks due to noise, since noise in these
%images is of the zero-mean type.
%find local maxima
fociCandidates = imregionalmax(IM,26);
%Find the mean of a local volume that will capture an entire point source.
xy = round(4*sigmaXY);
z = round(4*sigmaZ);
IMPad = padarray(IM,[xy xy z],'symmetric');
fociCandidatesPad = padarray(fociCandidates,[xy xy z]);
index = find(fociCandidatesPad);
clear fociCandidatesPad;
s=size(IMPad);
IMMeanIntensity = zeros(s);
if iscolumn(index)
    index = index';
end
for i = index
    [i2,j2,k2] = ind2sub(s,i);
    localIntensityVolume=IMPad(i2-xy:i2+xy,j2-xy:j2+xy,k2-z:k2+z);
    IMMeanIntensity(i) = mean(mean(mean(localIntensityVolume)));
end
clear IMPad; %clear up some memory
IMMeanIntensity = IMMeanIntensity(xy+1:end-xy,xy+1:end-xy,z+1:end-z);
%The final test statistic is the product of test statistic 1 and 2
spotStat = IMMeanIntensity(i).*curvature;
%----- Find a threshold that separates signal from noise -----
index = find(fociCandidates);
clear fociCandidates
if iscolumn(index)
    index = index';
end
spotStat = spotStat(index);
%find a good guess for the threshold using the triangle threshold
threshold = trithresh(spotStat);

if isrow(spotStat)
    spotStat = spotStat';
end
spotStat = [spotStat,index];
spotStat = sortRows(spotStat,1);
spotNum = size(spotStat,1);
for i = 1:spotNum
    if spotStat(i,1) > threshold
        ind = i;
        break;
    end
end

%search around the threshold space centered at the intial guess and find
%the threshold that minimizes the rate of change of foci count

s = size(IM);
foci = zeros(s);
for i = ind:spotNum
    foci(spotStat(i,2)) = 1;
end
%----- Create the final image with bonafide spots and other aesthetic images -----

%sum projection of foci
sumProj = sum(foci,3);

%max project the stamp
stampProj = zeros(size(sumProj));
%gaussian stamp (approx. the PSF) on the sum projection of foci
%3D Gaussian Filter\Stamp\PSF approximation
mu = [0,0,0]; %zero mean gaussian
SIGMA = [sigmaXY,0,0;0,sigmaXY,0;0,0,sigmaZ];
h = zeros(2*xy+1,2*xy+1,2*z+1);
for i=1:2*xy+1
    for j=1:2*xy+1
        for k=1:2*z+1
            h(i,j,k) = mvnpdf([i-1-xy,j-1-xy,k-1-z],mu,SIGMA); %the 3D filter
        end
    end
end
index = find(sumProj);
if iscolumn(index)
    index = index';
end
s = size(sumProj);
for i=index
    [i2,j2] = ind2sub(s,i);
 stampProj(i2-xy:i2+xy,j2-xy:j2+xy,k2-z:k2+z) = sumProj(i)*h;
end

for i = 1:lenz
    imwrite(uint16(IM3(:,:,i)),'fishtest6.tif','tif','WriteMode','append','Compression','none');
end

end

function [curvature] = mySobelHessianCurvature(I)
%This function uses the Sobel filter to approximate the derivatives in a
%gradient. Since the sobel filter is seperable the it can also be
%conveniently extended to find the Hessian.
%The image I is 3D and has coordinates (y,x,z).
h1 = [0.25 0.5 0.25];
h2 = [-1 0 1];
[~,sx,sz] = size(I);
tempI1 = zeros(size(I));
tempI2 = zeros(size(I));
%The Sobel separted filters to find the Fx
Fx = zeros(size(I));
for i=1:sx
    tempI1(:,i,:) = imfilter(I(:,i,:),h1,'symmetric'); %z
end
for i=1:sz
    tempI2(:,:,i) = imfilter(tempI1(:,:,i),h1,'symmetric'); %y
end
for i=1:sz
    Fx(:,:,i) = imfilter(tempI2(:,:,i),h2','symmetric'); %x
end
%The Sobel separted filters to find the Fy
Fy = zeros(size(I));
for i=1:sx
    tempI1(:,i,:) = imfilter(I(:,i,:),h1,'symmetric'); %z
end
for i=1:sz
    tempI2(:,:,i) = imfilter(tempI1(:,:,i),h2,'symmetric'); %y
end
for i=1:sz
    Fy(:,:,i) = imfilter(tempI2(:,:,i),h1','symmetric'); %x
end
%The Sobel separted filters to find the Fz
Fz = zeros(size(I));
for i=1:sx
    tempI1(:,i,:) = imfilter(I(:,i,:),h2,'symmetric'); %z
end
for i=1:sz
    tempI2(:,:,i) = imfilter(tempI1(:,:,i),h1,'symmetric'); %y
end
for i=1:sz
    Fz(:,:,i) = imfilter(tempI2(:,:,i),h1','symmetric'); %x
end
%The Sobel separted filters to find the Fxx
Fxx = zeros(size(I));
for i=1:sx
    tempI1(:,i,:) = imfilter(Fx(:,i,:),h1,'symmetric'); %z
end
for i=1:sz
    tempI2(:,:,i) = imfilter(tempI1(:,:,i),h1,'symmetric'); %y
end
for i=1:sz
    Fxx(:,:,i) = imfilter(tempI2(:,:,i),h2','symmetric'); %x
end
%The Sobel separted filters to find the Fxy
Fxy = zeros(size(I));
for i=1:sx
    tempI1(:,i,:) = imfilter(Fx(:,i,:),h1,'symmetric'); %z
end
for i=1:sz
    tempI2(:,:,i) = imfilter(tempI1(:,:,i),h1,'symmetric'); %y
end
for i=1:sz
    Fxy(:,:,i) = imfilter(tempI2(:,:,i),h2','symmetric'); %x
end
%The Sobel separted filters to find the Fxz
Fxz = zeros(size(I));
for i=1:sx
    tempI1(:,i,:) = imfilter(Fx(:,i,:),h1,'symmetric'); %z
end
for i=1:sz
    tempI2(:,:,i) = imfilter(tempI1(:,:,i),h1,'symmetric'); %y
end
for i=1:sz
    Fxz(:,:,i) = imfilter(tempI2(:,:,i),h2','symmetric'); %x
end
%The Sobel separted filters to find the Fyy
Fyy = zeros(size(I));
for i=1:sx
    tempI1(:,i,:) = imfilter(Fy(:,i,:),h1,'symmetric'); %z
end
for i=1:sz
    tempI2(:,:,i) = imfilter(tempI1(:,:,i),h2,'symmetric'); %y
end
for i=1:sz
    Fyy(:,:,i) = imfilter(tempI2(:,:,i),h1','symmetric'); %x
end
%The Sobel separted filters to find the Fyz
Fyz = zeros(size(I));
for i=1:sx
    tempI1(:,i,:) = imfilter(Fy(:,i,:),h1,'symmetric'); %z
end
for i=1:sz
    tempI2(:,:,i) = imfilter(tempI1(:,:,i),h2,'symmetric'); %y
end
for i=1:sz
    Fyz(:,:,i) = imfilter(tempI2(:,:,i),h1','symmetric'); %x
end
%The Sobel separted filters to find the Fzz
Fzz = zeros(size(I));
for i=1:sx
    tempI1(:,i,:) = imfilter(Fz(:,i,:),h2,'symmetric'); %z
end
for i=1:sz
    tempI2(:,:,i) = imfilter(tempI1(:,:,i),h1,'symmetric'); %y
end
for i=1:sz
    Fzz(:,:,i) = imfilter(tempI2(:,:,i),h1','symmetric'); %x
end
%Find the curvature matrix
clear tempI1;
clear tempI2;
curvature = zeros(size(I));
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
    curvature(i) = det(myHessian);
end
end

function [trithresh]=trithresh(A)
%----- Approximate histogram as triangle -----
%Create the histogram
[n,xout]=hist(A,100);
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
trithresh=xout(ind2);
end