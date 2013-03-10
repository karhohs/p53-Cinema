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
p.addParamValue('binning',1,@(x)any(bsxfun(@eq,x,[1,2])));
p.addParamValue('fluorchan','YFP',@(x)ischar(x));
p.parse(stackpath, varargin{:});
p2.met = p.Results.method;
p2.mag = p.Results.magnification;
p2.bin = p.Results.binning;

disp(['working in ', stackpath]); %Sanity Check
dirCon_stack = dir(stackpath);
stacknames=importStackNames(dirCon_stack,p.Results.fluorchan);

tempfoldername=regexp(stackpath,'(?<=\\)[\w ]*','match'); %Prepare to create a new folder to place background subtracted stacks
tempfoldername=[tempfoldername{end},'_bkgd'];
bkgdstackpath=[stackpath,'\..\',tempfoldername];
if ~exist(bkgdstackpath,'dir')
    mkdir(bkgdstackpath);
else
    disp('Note: Background substraction destination folder already exists.')
end
for i=1:length(stacknames)
    disp(['Subtracting background from ',stacknames{i}])
    info = imfinfo([stackpath,'\',stacknames{i}],'tif');
    t = Tiff([stackpath,'\',stacknames{i}],'r');
    Name = regexprep(stacknames{i},'(?<=_t)(\w*)(?=\.)','$1_bkgd');
    Name = fullfile(bkgdstackpath,Name);
    if length(info) > 1
        for k=1:length(info)-1
            IM = double(t.read);
            IM=bkgdmethods(IM,p2);
            IM = uint16(IM);
            imwrite(IM,Name,'tif','WriteMode','append','Compression','none');
            t.nextDirectory;
        end
        %one last time without t.nextDirectory
        IM = double(t.read);
        IM = bkgdmethods(IM,p2);
        IM = uint16(IM);
        imwrite(IM,Name,'tif','WriteMode','append','Compression','none');
    else
        IM = double(t.read);
        IM = bkgdmethods(IM,p2);
        IM = uint16(IM);
        imwrite(IM,Name,'tif','WriteMode','append','Compression','none');
    end
    %add image description from old stack to new stack
    t.setDirectory(1)
    metadata = t.getTag('ImageDescription');
    t.close;
    t = Tiff(Name,'r+');
    t.setTag('ImageDescription',metadata);
    t.rewriteDirectory;
    t.close;
end

end

function [S]=bkgdmethods(S,p)
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
switch lower(p.met)
    %Depending on the se2Size this method is similar to the 'uri' method.
    %It is the most conservative of all the methods. It is roughly 50 times
    %slower than the 'uri' method.
    case 'jared'
        switch p.mag
            %The values below have been heuristically chosen.
            case 10
                se1Size = 10;
                se2Size = 50;
            case 20
                se1Size = 20;
                se2Size = 100;
            case 40
                se1Size = 40;           % intracellular features in pixels, such as the nucleolus
                se2Size = 200;          % cells/cell clumps in pixels '' '' '
            case 60
                se1Size = 60;
                se2Size = 300;
            case 100
                se1Size = 100;
                se2Size = 500;
            otherwise
                error('bkgd:Whatever','How did you get here?');
        end
        switch p.bin
            case 1
                
            case 2
                se1Size = se1Size/2;
                se2Size = se2Size/2;
            otherwise
                error('bkgd:Whatever','How did you get here?');
        end
        resizeMultiplier = 1/2; % Downsampling scale factor makes image processing go faster and smooths image
        se1 = strel('disk', round(se1Size*resizeMultiplier));  %Structing elements are necessary for using MATLABS image processing functions
        se2 = strel('disk',round(se2Size*resizeMultiplier));
        
        origSize  = size(S);
        % Rescale image and compute background using closing/opening.
        I    = imresize(S, resizeMultiplier);
        % Pad image with a reflection so that borders don't introduce artifacts
        pad   = round(se2Size*resizeMultiplier);
        I    = padarray(I, [pad,pad], 'symmetric', 'both');
        % Perform opening/closing to get background
        I     = imclose(I, se1);   % ignore small low-intensity features (inside cells)
        I     = imopen(I, se2);     % ignore large high-intensity features (that are cells)
        % Remove padding and resize
        I     = floor(imresize(I(pad+1:end-pad, pad+1:end-pad), origSize));
        % Subtract background!
        S = S - I;
        S(S<0)=0;
    case 'alex'
        %This algorithm is more aggressive than the 'uri' algorithm. It
        %also takes approximately 100 times longer to run then the 'uri'
        %algorithm. I do not think this works well at large grid
        %sizes/high-magnification
        
        switch p.mag
            %The values below have been heuristically chosen.
            case 10
                sel=50;
            case 20
                sel=100;
            case 40
                sel=200;
            case 60
                sel=300;
            case 100
                sel=500;
            otherwise
                error('bkgd:Whatever','How did you get here?');
        end
        switch p.bin
            case 1
                
            case 2
                sel = sel/2;
            otherwise
                error('bkgd:Whatever','How did you get here?');
        end
        
        bkgd=gridscan(S,@gaussianthreshold,round(sel));
        S=S-bkgd;
        S(S<0)=0;
        
    case 'uri'
        %This algorithm seems to work well when the sel value is about the
        %size of a large nucleus. Larger senescent nuclei might prove
        %troublesome if the sel is based on large nuclei at the start of
        %the movie, which would be considerably smaller. The larger the sel
        %the more conservative the estimate of the background. Between the
        %'jared', 'alex', and 'uri' methods the 'uri' method is the fastest
        %by at least 2 orders of magnitude.
        
        %I do not think this works well at large-grid-sizes/high-magnification
        switch p.mag
            %The values below have been heuristically chosen.
            case 10
                sel=50;
            case 20
                sel=100;
            case 40
                sel=200;
            case 60
                sel=300;
            case 100
                sel=500;
            otherwise
                error('bkgd:Whatever','How did you get here?');
        end
        switch p.bin
            case 1
                
            case 2
                sel = sel/2;
            otherwise
                error('bkgd:Whatever','How did you get here?');
        end
        bkgd=gridscan(S,@rankfilter,sel);
        S=S-bkgd;
        S(S<0)=0;
    case 'tonia'
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

