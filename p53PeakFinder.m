function [pks,vly,pksPitch,vlyPitch,sus,susPitch] = p53PeakFinder(signal,time,samplingFreq)
% Input:
% signal: presumabley a vector of time varying fluorescent protein data
% time: the times at which measurements were taken
% samplingFreq:
%
% Output:
% pks: the time points where peaks in the data exist
% vly: the time points where valleys in the data exist
% pksPitch: the pitch for each peak. A measure of how dense peaks are
% surrounding this peak.
% vlyPitch: the pitch for each valley. A measure of how dense valleys are
% surrounding this valley. Note: valley data should closely resemble the
% peak data.
% sus: the time points where a switch in protein levels occurs. Sus is
% short for sustained, because of the desire to find the time point where
% protein levels switch from low expression to high expression.
% susPitch: This measure should give a sense of how long sustained
% expression lasts as a measure of time. Sustained expression can be
% thought of as very long plateaued pulses.
%
% Description:
%
%
% Other Notes:
% waveinfo('mexh')
% To locate both peaks, valleys, and their pitch use mexican hat wavelet
% waveinfo('gaus') derivativeOfGaussian degree 1
% To locate singularities and switches use the derivativeOfGaussian degree 1 wavelet.

[s2,t] = scrubData(signal,time);
%Baseline removal: Remove the baseline using the discrete wavelet transform
%try imopen(x,strel('disk',w))

%Find the de-scaled continous wavelet transform. De-scaled refers to
%undoing the scaling done to ensure energy preservation.
wavMexh = cwtftNonuniformScalesMexh(s2);
wavDog1 = cwtftNonuniformScalesDog1(s2);

%remove padding from signal and wavelet transforms to only analyze the true signal.
temp = length(s2)-length(t);
if mod(temp,2)
    %is odd
    temp=temp-1;
    s = s2(temp/2+1:end-temp/2-1);
    wavMexh.cfs = wavMexh.cfs(:,temp/2+1:end-temp/2-1);
    wavDog1.cfs = wavDog1.cfs(:,temp/2+1:end-temp/2-1);
else
    %is even
    s = s2(temp/2+1:end-temp/2);
    wavMexh.cfs = wavMexh.cfs(:,temp/2+1:end-temp/2);
    wavDog1.cfs = wavDog1.cfs(:,temp/2+1:end-temp/2);
end

%Find the ridgemap for peak detection using the continuous wavelet
%transform.
x = findRidgeMap(wavMexh.cfs)
end

function [s,t] = scrubData(signal,time)
%Interpolate the signal to increase the density of points to increase
%sensitivity to high frequency content, to ensure data points are evenly
%spaced in time, and to fill in time points to a power of 2 to make use of
%the FFT for wavelet transformation.
t_diff = diff(time); %find time interval
temp = median(t_diff)/3; % This equation is a sort of "secret sauce". Practically speaking, images will be taken at fixed intervals, but occassionally there are outlier intervals due to things such as irradiating cells. Therefore the median should filter out outliers and hone in on the typical behavior. Finally, dividing by 3 will increase the density of data points by 3 using interpolation. This helps identify high frequency content. 3 was chosen empirically.
t = (time(1):temp:time(end)); %time points are equally spaced
s = spline(time,signal,t); %Interpolate with splines
% Find the next power of two that greater than 150% the current.
fftLen= 2^(ceil(log(length(s)*1.5)/log(2)));
%Extend signal to fill the FFT length, with a slight variation for filling
%even or odd number of padding
temp = fftLen-length(s);
if mod(temp,2)
    %is odd
    s = padarray(s,[0 (temp-1)/2],'symmetric');
    s(end+1) = s(end);
else
    %is even
    s = padarray(s,[0 temp/2],'symmetric');
end
%The wavelet transform is sensitive to edge effects. Since the signal
%(almost always) begins and ends abruptly, i.e. not equal to zero, the cone
%of influence from the edges will bleed into signal at high scales. To
%reduce this unwanted edge effect a tukey window is applied to the input
%signal.
tukey33=tukeywin(length(s),0.33)';
s = s.*tukey33;
end

function [out] = findRidgeMap(in)
%Input:
%in: a heatmap of the continous wavelet transform coefficients
%
%Output:
%out: a cell array with ridges within each cell. Each ridge is a 2xn vector
%where n is the length the the ridge. The first row of the vector contains
%the scale component and the second row contains the time component.

%-- Step 1: Perform Continuous Wavelet Transform --
%note: input signal 'y' must be evenly spaced
wavelet_scales = [1,3,5,7,(9:floor(length(y)/10):length(y))];
%wavelet_xfrm_coefs = cwt(y,wavelet_scales,'haar');
wavelet_xfrm_coefs = cwt(y,wavelet_scales,'mexh');
wavelet_peaks = cell(size(wavelet_xfrm_coefs,1),1);
for i=1:size(wavelet_xfrm_coefs,1)
    wavelet_peaks{i} = first_pass_peak_detection(...
        wavelet_xfrm_coefs(i,:),wavelet_scales(i)*2+1); %The window size was heursitically chosen to be (2*current_scale+1)
end
%-- Identify the Ridges --
ridge_map = zeros(size(wavelet_xfrm_coefs));
wavelet_peaks_alias = cell(size(wavelet_peaks));
gap_limit = 2;

for i=1:length(wavelet_peaks_alias)
    wavelet_peaks_alias{i} = zeros(size(wavelet_peaks{i}));
end
%The ridges are first seeded with the peaks in the highest scale
for i=1:length(wavelet_peaks{end})
    ridge_map(end,wavelet_peaks{end}(i)) = i;
    wavelet_peaks_alias{end}(i) = i;
end
ridge_counter = length(wavelet_peaks{end}); %keeps track of the total number of ridges
for i=size(wavelet_xfrm_coefs,1):-1:(gap_limit+1)
    
    for j=1:length(wavelet_peaks{i})
        for h=1:gap_limit
            %Search for peaks within the window size for scale i.
            low_bnd = wavelet_peaks{i}(j) - wavelet_scales(i-h);
            up_bnd = wavelet_peaks{i}(j) + wavelet_scales(i-h);
            low_set = wavelet_peaks{i-h}>low_bnd;
            up_set = wavelet_peaks{i-h}<up_bnd;
            %If a peak is found add it to the growing ridge
            if any(low_set.*up_set)
                ridge_set = low_set.*up_set;
                if sum(ridge_set)==1
                    for k=1:length(ridge_set)
                        if ridge_set(k) && (wavelet_peaks_alias{i-h}(k)==0 && ...
                                sum(wavelet_peaks_alias{i-h}==wavelet_peaks_alias{i}(j))==0)
                            ridge_map(i-h,wavelet_peaks{i-h}(k)) = wavelet_peaks_alias{i}(j);
                            wavelet_peaks_alias{i-h}(k) = wavelet_peaks_alias{i}(j);
                        end
                    end
                else
                    my_min = inf;
                    for k=1:length(ridge_set)
                        if ridge_set(k) && wavelet_peaks_alias{i-h}(k)==0
                            temp_ind=wavelet_peaks_alias{i-(h-1)}==wavelet_peaks_alias{i}(j);
                            penultimate_ridge_position=wavelet_peaks{i-(h-1)}(temp_ind);
                            temp_min = abs(wavelet_peaks{i-h}(k)-penultimate_ridge_position);
                            if temp_min<my_min
                                k_min = k;
                                my_min = temp_min;
                            end
                        end
                    end
                    if sum(wavelet_peaks_alias{i-h}==wavelet_peaks_alias{i}(j))==0
                        ridge_map(i-h,wavelet_peaks{i-h}(k_min)) = wavelet_peaks_alias{i}(j);
                        wavelet_peaks_alias{i-h}(k_min) = wavelet_peaks_alias{i}(j);
                    end
                end
            end
        end
    end
    %Start new ridges for the peaks that were not assigned to ridges in the
    %scale below.
    new_ridge_set = wavelet_peaks_alias{i-1}==0;
    for j=1:length(wavelet_peaks_alias{i-1})
        if new_ridge_set(j)
            ridge_counter = ridge_counter + 1;
            wavelet_peaks_alias{i-1}(j) = ridge_counter;
            ridge_map(i-1,wavelet_peaks{i-1}(j)) = ridge_counter;
        end
    end
end

%Find the peaks of every ridge
ridge_peaks=zeros(ridge_counter,2);
%weighted_ridge_peaks=zeros(ridge_counter,2);
for i=1:ridge_counter
    temp_map=ridge_map==i;
    temp_map2=wavelet_xfrm_coefs-min(min(wavelet_xfrm_coefs));
    temp_map2(~temp_map)=0;
    [temp_max_vector,ind_vector]=max(temp_map2);
    [~,ind]=max(temp_max_vector);
    ridge_peaks(i,:)=[ind,ind_vector(ind)];
    
    %    temp_map=ridge_map==i;
    %    temp_map2=weighted_wavelet_xfrm_coefs;
    %    temp_map2(~temp_map)=0;
    %    [temp_max_vector,ind_vector]=max(temp_map2);
    %    [~,ind]=max(temp_max_vector);
    %    weighted_ridge_peaks(i,:)=[ind,ind_vector(ind)];
end

%The peaks of every ridge represent a candidate peak from the original
%waveform. The ridge peak contains both positional and scale information.
%Selective criteria based upon the scale can be used to sift through noise
%and get a sense for the breadth of each peak.
% peaks=ridge_peaks(:,1); %3 is chosen heursitically as the noise floor.
% peaks=peaks(ridge_peaks(:,2)>3);

%Choose peaks that are on a scale that best picks up p53 pulses
peaks_p53_pulses=find(ridge_map(5,:)~=0);
peaks=peaks_p53_pulses;
end

function [peak_index]=first_pass_peak_detection(x,window_size)
%There are many ways to find peaks. The wavelet method is powerful at
%detecting peaks but is actually dependent on simpler peak detection
%methods. A simple approach that is itself effective is to identify
%candidate peaks and then weed out unwanted peaks through thresholds or
%selective criteria.
%Step 1: Indentify candidate peaks.
%Method 1: Use the watershed algorithm. If the signal is sufficiently
%smooth this approach is awesome.
a=2;
if a == 1
    peak_index = watershed(x);
    peak_index = find(~peak_index);
else
    %Method 2: Scan a signal with a window. Find the max. If the max is greater
    %than the left and right endpoints of the window it is a peak candidate.
    %The window is centered at each point of a waveform, so a wider peak
    %candidate will recieve more votes in a sense. If a window has more than
    %one point with the max value then the left most index is used. A downside
    %to this method is that it is sensitive to the size of the window. However,
    %this seeming limitation is actually put to use in the wavelet method by
    %scaling the window along with the wavelet.
    peak_index = window_peak_finder(window_size,x);
end

%Weed out unqualified peaks and peaks of questionable nature.

end

function [out]=window_peak_finder(winw_size,x)
length_x = length(x);
if size(x,2) == 1
    x=x';
end

pad=(winw_size-1)/2;
[x_padded,window_padded] = deal(zeros(1,length_x+winw_size-1));
window_padded(1:winw_size) = ones(1,winw_size);
x_padded(pad+1:end-pad) = x;
out=zeros(size(x));
for i=1:length_x
    seg = x_padded(logical(window_padded));
    [~,ind] = max(seg);
    if seg(ind) > seg(1) && seg(ind) > seg(end)
        out(i) = i + ind - pad - 1;
    end
    window_padded = circshift(window_padded,[0, 1]);
end
out(out==0) = [];
peak_candidates = unique(out);
peak_elected = peak_candidates; %This vector will be trimmed below
%Tally votes; if a peak candidate has less than or equal to pad votes they
%are disqualified
for i=length(peak_candidates):-1:1
    temp = out==peak_candidates(i);
    if sum(temp)<=pad/2
        peak_elected(i) = [];
    end
end
%If a peak has a negative wavelet coefficient value it is disqualified.
for i=length(peak_elected):-1:1
    if x(peak_elected(i))<0
        peak_elected(i) = [];
    end
end

out = peak_elected;
end

function [my_fig]=plot_peaks(peak_index,x)
my_fig = figure;
plot(x)
hold on
x2=x;
temp=(1:length(x));
temp(peak_index)=[];
x2(temp)=NaN;
plot(x2,'x','color','r')
hold off
end

function []=rapid_report(wavelet_peaks,wavelet_xfrm_coefs,ridge_map,x,y,peaks_p53_pulses)
temp = clock;
time_stamp = regexprep(num2str(temp(1:5)),'\s','');
filename = ['wavelet_analysis_' time_stamp];

my_fig = figure;
plot(x,y)
title('Original Waveform')
print(my_fig,[filename '.ps'],'-dpsc2','-painters');
close(my_fig)

my_fig = plot_peaks(wavelet_peaks{1},wavelet_xfrm_coefs(1,:));
my_fig_name = 'Wavelet Coefficients, Scale 1';
title(my_fig_name)
print(my_fig,[filename '.ps'],'-dpsc2','-painters','-append');
close(my_fig)
for i=2:length(wavelet_peaks)
    my_fig = plot_peaks(wavelet_peaks{i},wavelet_xfrm_coefs(i,:));
    temp = num2str(2*i-1);
    my_fig_name = ['Wavelet Coefficients, Scale ' temp];
    title(my_fig_name)
    print(my_fig,[filename '.ps'],'-dpsc2','-painters','-append');
    close(my_fig)
end

my_fig = figure;
imagesc(ridge_map)
colormap('jet')
title('Ridge Map')
print(my_fig,[filename '.ps'],'-dpsc2','-painters','-append');
close(my_fig)

my_fig = figure;
plot(x,y,'Color','black')
hold on
temp_y=NaN(size(y));
temp_peaks=y(peaks_p53_pulses);
temp_y(peaks_p53_pulses)=temp_peaks;
plot(x,temp_y,'d','color','red','MarkerSize',6,'MarkerFaceColor','red')
title('Original Waveform with Peaks')
hold off
print(my_fig,[filename '.ps'],'-dpsc2','-painters','-append');
close(my_fig)

ps2pdf('psfile', [filename '.ps'], 'pdffile', [filename '.pdf'], 'gspapersize', 'a4', 'deletepsfile', 1);
end

function [out] = cwtftNonuniformScalesMexh(in)
%Input:
%in: the 1D time-domain signal
%
%Output:
%out: a struct with the wavelet transformation.
%out.cfs = wavelet coefficients in a matrix
%out.scl = the scales of the wavelet
%out.phz = the pseudofrequencies of the wavelet scales
%
%Description:
%The cwtft function will only calculate coefficients
%for evenly spaced scales. However, this is an inconvenience when trying to
%look at trends that occur at different scales. To overcome this difficulty
%this function calls the cwtft function iteratively and then assembles the
%data into a struct that contains the coefficients and scales.

%Choosing the right scales to investigate can be a challenge, because it
%can feel subjective. One way is to include every integer scale up to the
%length of the signal, but this is probably too much information. Another
%is to choose scales on an exponential/log scale, but this might gloss over
%some important details. I will try some hybrid between the two. Filling in
%between an exponential scale with uniform spacing. Hopefully this
%comprimise will deliver detail across several orders of magnitude.

%It was found to be that wavelet coefficients are no longer useful once
%the wavelet support is approx. half the length of the signal. The Mexican
%Hat wavelet has a support of 8 (or is it 11? check waveinfo('mexh')) at
%scale 1. Therefore, in order to find out where the "half support is"...
if length(in)<=32
    warning('wavylsis:tooshort', 'The length of the input signal may be too short for wavelet analysis');
end
hs = length(in)/16;
%Now implement a scale scheme
pow102 = ceil(log(hs/10)/log(2));
if pow102<2
    wavelet_scales{1} = (1:10);
    wavelet_scales{pow102}(wavelet_scales{pow102}>hs) = [];
    temp = cwtft(in,'scales',wavelet_scales{1},'wavelet','mexh');
    out.cfs = temp.cfs;
    out.scl = wavelet_scales{1};
    out.phz = wavelet_scales{1}*0.25; %0.25 is the pseudofrequency of the mexh for scale 1
else
    wavelet_scales = cell(1,pow102);
    for i=1:pow102
        wavelet_scales{i} = (1:10)*2^(i-1)+10*(2^(i-1)-1);
    end
    wavelet_scales{pow102}(wavelet_scales{pow102}>hs) = [];
    temp = cwtft(in,'scales',wavelet_scales{1},'wavelet','mexh');
    out.cfs = temp.cfs;
    for i=2:pow102
        temp = cwtft(in,'scales',wavelet_scales{i},'wavelet','mexh');
        temp = temp.cfs;
        out.cfs = [out.cfs;temp];
    end
    
    out.scl = cell2mat(wavelet_scales);
    out.phz = out.scl*0.25; %0.25 is the pseudofrequency of the mexh for scale 1
end
out.cfs = real(out.cfs);
for i = 1:length(out.scl)
    out.cfs(i,:) = out.cfs(i,:)/sqrt(out.scl(i));
end
end

function [out] = cwtftNonuniformScalesDog1(in)
%Input:
%in: the 1D time-domain signal
%
%Output:
%out: a struct with the wavelet transformation.
%out.cfs = wavelet coefficients in a matrix
%out.scl = the scales of the wavelet
%out.hz = the pseudofrequencies of the wavelet scales
if length(in)<=32
    warning('wavylsis:tooshort', 'The length of the input signal may be too short for wavelet analysis');
end
hs = length(in)/16;
%Now implement a scale scheme
pow102 = ceil(log(hs/10)/log(2));
if pow102<2
    wavelet_scales{1} = (1:10);
    wavelet_scales{pow102}(wavelet_scales{pow102}>hs) = [];
    temp = cwtft(in,'scales',wavelet_scales{1},'wavelet',{'dog',1});
    out.cfs = temp.cfs;
    out.scl = wavelet_scales{1};
    out.phz = wavelet_scales{1}*0.2; %0.2 is the pseudofrequency of the dog1 for scale 1
else
    wavelet_scales = cell(1,pow102);
    for i=1:pow102
        wavelet_scales{i} = (1:10)*2^(i-1)+10*(2^(i-1)-1);
    end
    wavelet_scales{pow102}(wavelet_scales{pow102}>hs) = [];
    temp = cwtft(in,'scales',wavelet_scales{1},'wavelet',{'dog',1});
    out.cfs = temp.cfs;
    for i=2:pow102
        temp = cwtft(in,'scales',wavelet_scales{i},'wavelet',{'dog',1});
        temp = temp.cfs;
        out.cfs = [out.cfs;temp];
    end
    out.scl = cell2mat(wavelet_scales);
    out.phz = out.scl*0.2; %0.2 is the pseudofrequency of the dog1 for scale 1
end
out.cfs = real(out.cfs);
for i = 1:length(out.scl)
    out.cfs(i,:) = out.cfs(i,:)/sqrt(out.scl(i));
end
end