function []=backgroundRemoval(stackpath,varargin)
% [] = backgroundRemoval()
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
% 
p = inputParser;
p.addRequired('stackpath', @(x)ischar(x));
p.addParamValue('method','Jared',@(x)ischar(x));
p.addParamValue('magnification',20,@(x)any(bsxfun(@eq,x,[10,20,40,60,100])));
p.addParamValue('binning',1,@(x)any(bsxfun(@eq,x,[1,2,4,8,16])));
p.parse(stackpath, varargin{:});
stacknames=importStackNames(stackpath);
tempfoldername=regexp(stackpath,'(?<=\\)[\w ]*','match'); %Prepare to create a new folder to place background subtracted stacks
tempfoldername=[tempfoldername{end},'_bkgd'];
bkgdstackpath=[stackpath,'\..\',tempfoldername];
if ~exist(bkgdstackpath,'dir')
    mkdir(bkgdstackpath);
    for i=1:length(stacknames)
        disp(['Subtracting background from ',stacknames{i}])
        S=readStack([stackpath,'\',stacknames{i}]);
        S=bkgdmethods(S,'jared');
        newStack(S,stacknames{i},bkgdstackpath);
    end
else
    disp('Note: Background subtraction has occurred previously.')
end
end

function [S]=bkgdmethods(S,method)
%S is a stack, method chooses between several different background
%subtraction methods
%Jared's method is a morphological approach to find the background. An
%opening and closing procedure is performed. The key is choosing
%appropriate structing element sizes. By default they are (16,close) and
%(320,open). These are chosen to be approximately the size of intracellular
%noise and clumps of cells respectively. Using (96,open) gives results
%closer to the other methods using a 32 pixel grid.
%
%Alex's method uses Gaussian fitting to identify local thresholds. The
%assumption here is that the background is well approximated by a Gaussian
%distribution and the signal of interest is significantly higher than this
%distribution, i.e. 2 to 3 sigma away. A local area is defined roughly 1 to
%2 times the size of a nucleus (in this case 96 pixels is chosen).
%
%Uri's method is based upon background subtraction perfromed by Sigal et
%al., Nature Methods 2006. The image is broken up into a grid and from each
%grid space the 10th percentile pixel value is chosen to as the background.
%The size of the grid must be roughly 1 to 2 times the size of a nucleus if
%measuring a nuclear localized protein (in this case 96 pixels is chosen).
switch lower(method)
    %Depending on the se2Size this method is similar to the 'uri' method.
    %It is the most conservative of all the methods. It is roughly 50 times
    %slower than the 'uri' method.
    case 'jared'
        resizeMultiplier = 1/2; % Downsampling scale factor makes image processing go faster and smooths image
        %se1Size = 20;           % Upper limit of intracellular features in pixels for a 672x512pixel image from a 40x objective inside MCF7 cells
        %se2Size = 200;          % Upper limit of cells/cell clumps in pixels '' '' ''
        %se1Size = 4;           % Upper limit of intracellular features in pixels for a 672x512pixel image from a 20x objective inside MCF7 cells
        %se2Size = 30;          % Upper limit of cells/cell clumps in pixels '' '' ''      
        se1Size = 4;           % Upper limit of intracellular features in pixels for a 1344x1024pixel image from a 60x objective inside FISH cells
        se2Size = 40;          % Upper limit of cells/cell clumps in pixels '' '' ''            
        se1 = strel('disk', se1Size*resizeMultiplier);  %Structing elements are necessary for using MATLABS image processing functions
        se2 = strel('disk',se2Size*resizeMultiplier);
        if iscell(S)
            origSize  = size(S{1});
            for k=1:length(S)
                % Rescale image and compute background using closing/opening.
                
                I    = imresize(S{k}, resizeMultiplier);
                pad   = round(se2Size*resizeMultiplier);
                % Pad image with a reflection so that borders don't introduce artifacts
                I    = padarray(I, [pad,pad], 'symmetric', 'both');
                % Perform opening/closing to get background
                I     = imclose(I, se1);   % ignore small low-intensity features (inside cells)
                I     = imopen(I, se2);     % ignore large high-intensity features (that are cells)
                % Remove padding and resize
                I     = floor(imresize(I(pad+1:end-pad, pad+1:end-pad), origSize));
                % Subtract background!
                S{k} = S{k} - I;
                S{k}(S{k}<0)=0;
                %NOTE: b/c I is an unsigned integer array, negative numbers
                %created by subtraction all become zero.
                
            end
        else
            origSize  = size(S);
            % Rescale image and compute background using closing/opening.
            I    = imresize(S, resizeMultiplier);
            pad   = round(se2Size*resizeMultiplier);
            I    = padarray(I, [pad,pad], 'symmetric', 'both');
            I     = imclose(I, se1);   % ignore small low-intensity features (inside cells)
            I     = imopen(I, se2);     % ignore large high-intensity features (that are cells)
            I     = floor(imresize(I(pad+1:end-pad, pad+1:end-pad), origSize));
            S = S - I;
            S(S<0)=0;
        end
    case 'alex'
        %This algorithm is more aggressive than the 'uri' algorithm. It
        %also takes approximately 100 times longer to run then the 'uri'
        %algorithm
        sel=100; %structuring element length. The 512x672 image will be broken into 32x32 squares
        if iscell(S)
            for i=1:length(S)
                bkgd=gridscan(S{i},@gaussianthreshold,sel);
                S{i}=S{i}-bkgd;
                S{i}(S{i}<0)=0;
            end
        else
            bkgd=gridscan(S,@gaussianthreshold,sel);
            S=S-bkgd;
            S(S<0)=0;
        end
    case 'uri'
        %This algorithm seems to work well when the sel value is about the
        %size of a large nucleus. Larger senescent nuclei might prove
        %troublesome if the sel is based on large nuclei at the start of
        %the movie, which would be considerably smaller. The larger the sel
        %the more conservative the estimate of the background. Between the
        %'jared', 'alex', and 'uri' methods the 'uri' method is the fastest
        %by at least 2 orders of magnitude.
        sel=300; %structuring element length. %100 for 20x, 200 for 40x, 300 for 60x
        if iscell(S)
            for i=1:length(S)
                bkgd=gridscan(S{i},@rankfilter,sel);
                S{i}=S{i}-bkgd;
                S{i}(S{i}<0)=0;
            end
        else
            bkgd=gridscan(S,@rankfilter,sel);
            S=S-bkgd;
            S(S<0)=0;
        end
    otherwise
        error('Unknown method of background subtraction specified')
end
end

function [bkgd]=gridscan(S,funfcn,sel)
%The pass over the inner grid
ho=size(S,1); %original height of the input image (usually 512)
wo=size(S,2); %original width of the input image (usually 672)

%ensure the image size is divisible by the grid size
xtrah=mod(ho,sel);
if xtrah
    h=ho+(sel-xtrah);
    S=padarray(S,[(sel-xtrah) 0],'symmetric','post');
else
    h=ho;
end
xtraw=mod(wo,sel);
if xtraw
    w=wo+(sel-xtraw);
    S=padarray(S,[0 (sel-xtraw)],'symmetric','post');
else
    w=wo;
end

gridh=h/sel;
gridw=w/sel;
bkgdgridcenter=zeros(gridh,gridw); %Holds the "center" background intensity value for each grid space

%Pass over the inner-grid
jarray=(sel:sel:h);
karray=(sel:sel:w);
for j=jarray
    for k=karray
        A=S(j-sel+1:j,k-sel+1:k);
        bkgdgridcenter((j)/sel,(k)/sel)=funfcn(A);
    end
end

%Create background image using interpolation with the bkgdgridcenter
bkgd=floor(imresize(bkgdgridcenter,sel));
%Remove padding if there was any
bkgd=bkgd(1:ho,1:wo);
end

function [B]=rankfilter(A)
B=reshape(A,[],1);
B=sort(B); %Sort values in the grid in ascending order. In other words rank the pixels by intensity.
%The percentile can be tweaked and adjusted to improve results
B=B(floor(length(B)*(0.1))); %Find the 10th percentile. The idea is that the 10th percentile is always representative of the background.
end

function [B]=gaussianthreshold(A)
%This function will attempt to fit up to three Gaussians. The mean of the
%first Guassian is taken to be the background threshold
%As of writing this code the documentation for the fitting toolbox is
%awful. I find trying to fit more than two guassians leads to funny fits,
%so I recommend not fitting more than two. Each Gaussian has three
%parameters: a=amplitude, b=mean, c=variance.
ftype = fittype('gauss1');
opts = fitoptions('gauss1');
A=reshape(A,[],1);
s=min(A);
r=max(A);
x = (s:floor((r-s)/40):r)';
h = hist(A,x)';
[FitModel, Goodness] = fit(x,h,ftype,opts);
if Goodness.rsquare>0.96
    B=floor(FitModel.b1);
else
    ftype = fittype('gauss2');
    opts = fitoptions('gauss2');
    [FitModel, ~] = fit(x,h,ftype,opts);
    B=floor(FitModel.b1);
    
end
end
