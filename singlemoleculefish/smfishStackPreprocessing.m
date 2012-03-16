%noise reduction by matched filtering with a guassian kernel
%create gaussian filter that approximates 3D PSF of the microscope.
%Distance units are in micrometers.
objective = 60; %as in 60x
NA = 1.4; %typical of 60x oil immersion objections
rindex = 1.51; %typical refractive index of oil
camerapixellength = 6.45; %Both cameras in the Lahav have pixel dimensions of 6.45 x 6.45 um.
zstepsize = 0.25; %User defined with z-stack is obtained
wavelength = .67; %Cy5 probe wavelength approximately 670 nanometers
sigmaXYos = 0.21*wavelength/NA; %lateral st. dev of the gaussian filter in object space
sigmaZos = 0.66*wavelength*rindex/(NA^2) ;%axial st. dev of the gaussian filter in object space
Pxy = camerapixellength/objective; %lateral pixel size
sigmaXY = sigmaXYos/Pxy; %lateral st. dev of gaussian filter in image space
sigmaZ = sigmaZos/zstepsize; %axial st. dev of gaussian filter in image space
xy = round(3*sigmaXY);
z = round(3*sigmaZ);
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
%apply filter to the image
path = 'M:\SI_201203B\deconvolution_test\orignal_data_ff_bkgd\20120312_pos1-1_w2Cy5_s1_t1_ff.TIF';
info = imfinfo(path,'tif');
t = Tiff(path,'r');
IM_temp = t.read;
[leny,lenx] = size(IM_temp);
lenz = length(info);
IM = zeros(leny,lenx,lenz);
if lenz > 1
    for k=1:length(info)-1
        IM(:,:,k) = uint16(t.read); %double(t.read);
        t.nextDirectory;
    end
    %one last time without t.nextDirectory
    IM (:,:,lenz) = uint16(t.read); %double(t.read);
end
%IM2 = imfilter(IM,h,'symmetric'); %Totally works. sweet!
%enhance spots. OMG Everything above is unnecessary :(
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
h=-h;
IM2 = imfilter(IM,h,'symmetric'); %Totally works. sweet!
for i = 1:lenz
imwrite(uint16(IM2(:,:,i)),'fishtest3.tif','tif','WriteMode','append','Compression','none');
end