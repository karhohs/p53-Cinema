function [pks,vly,sus,ssas,maps] = p53PeakFinder(signal,time,samplingFreq)
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
% ssas: the single-sided amplitude spectrum from the fourier transform of the signal
% maps: A collection of 2D matrices
% maps.mexh.ridg
% maps.mexh.cwtcfs
% maps.dog1.ridg
% maps.dog1.cwtcfs
%
% Description:
%
%
% Other Notes:
% waveinfo('mexh')
% To locate both peaks, valleys, and their pitch use mexican hat wavelet
% waveinfo('gaus') derivativeOfGaussian degree 1
% To locate singularities and switches use the derivativeOfGaussian degree 1 wavelet.

%Find the ssas
%ssas = findSSAS(signal,time,samplingFreq);

%Baseline removal: Remove the baseline using the discrete wavelet transform
%try imopen(x,strel('disk',w))
sbr = imopen(signal,strel('disk',5));

[s2,t] = scrubData(signal,time);

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
%transform. The signal is given to the findRidgeMap to use as a threshold.
[ridgpks,wavpeaks] = findRidgeMap(wavMexh,s);
if isempty(ridgpks.scl)
    disp('There are no peaks in this waveform')
    return
end
%opportunity here to prune ridges by length
[protopks,pkSclStats] = processRidgeMap(ridgpks,wavMexh.scl,wavMexh.phz,length(s));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Everything above here is pretty solid
%Repeat ridge mapping for valleys and switches
temp = wavMexh;
temp.cfs = -temp.cfs;
[ridgval,waveval] = findRidgeMap(temp,s);
[protoval,valSclStats] = processRidgeMap(ridgval,wavMexh.scl,wavMexh.phz,length(s));
beautifyRidgeMap(ridgpks.map,ridgval.map,wavMexh.cfs);

%switches
[ridgeSwitchPks,waveSwitchPks] = findRidgeMap(wavDog1,s);
[protoSwitchPks,pkSwitchSclStats] = processRidgeMap(ridgeSwitchPks,wavDog1.scl,wavDog1.phz,length(s));

%%%%%%%%%%%%%%
ptime = cell(1,length(protoSwitchPks));
[ptime{:}] = protoSwitchPks.time;
ptime = cell2mat(ptime);
psclind = cell(1,length(protoSwitchPks));
[psclind{:}] = protoSwitchPks.scaleindex;
psclind = cell2mat(psclind);
palltogether = [ptime',psclind'];
newmap = zeros(size(ridgpks.map));
for i=1:length(palltogether)
    j = palltogether(i,2);
    k = palltogether(i,1);
    newmap(j,k) = 1;
end
beautifyRidgeMap(ridgeSwitchPks.map,newmap,wavDog1.cfs);
for i=1:length(palltogether)
    j = palltogether(i,2);
    k = palltogether(i,1);
    newmap(j,k) = i;
end

figure
imagesc(newmap);
penguinjet;
%The peaks of every ridge represent a candidate peak from the original
%waveform. The ridge peak contains both positional and scale information.
%Selective criteria based upon the scale can be used to sift through noise
%and get a sense for the breadth of each peak.
pks = processPeaks(s,protopks,pkSclStats);
vly = processPeaks(-s,protoval,valSclStats);
%<DEBUG>
ptime = cell(1,length(protoval));
[ptime{:}] = protoval.time;
ptime = cell2mat(ptime);
pscl = cell(1,length(protoval));
[pscl{:}] = protoval.scale;
pscl = cell2mat(pscl);
pcfs = cell(1,length(protoval));
[pcfs{:}] = protoval.waveletcfs;
pcfs = cell2mat(pcfs);
plen = cell(1,length(protoval));
[plen{:}] = protoval.ridgelength;
plen = cell2mat(plen);
valltogether = [ptime',pscl',pcfs',plen'];
valltogether = sortrows(valltogether,1);

ptime = cell(1,length(protopks));
[ptime{:}] = protopks.time;
ptime = cell2mat(ptime);
pscl = cell(1,length(protopks));
[pscl{:}] = protopks.scale;
pscl = cell2mat(pscl);
psclind = cell(1,length(protopks));
[psclind{:}] = protopks.scaleindex;
psclind = cell2mat(psclind);
pcfs = cell(1,length(protopks));
[pcfs{:}] = protopks.waveletcfs;
pcfs = cell2mat(pcfs);
plen = cell(1,length(protopks));
[plen{:}] = protopks.ridgelength;
plen = cell2mat(plen);
ptype = cell(1,length(pks));
[ptype{:}] = pks.type;
ptype = cell2mat(ptype);
palltogether = [ptime',pscl',pcfs',plen',ptype',psclind'];
palltogether = sortrows(palltogether,1);
newmap = zeros(size(ridgpks.map));
for i=1:length(palltogether)
    j = palltogether(i,6);
    k = palltogether(i,1);
    newmap(j,k) = 1;
end
beautifyRidgeMap(ridgpks.map,newmap,wavMexh.cfs);
%</DEBUG>
rapid_report(waveval,-wavMexh.cfs,ridgval.map,t,s)
plotPeaksAndValleys(s,pks,vly)
end

function [s,t] = scrubData(signal,time)
%Interpolate the signal to increase the density of points to increase
%sensitivity to high frequency content, to ensure data points are evenly
%spaced in time, and to fill in time points to a power of 2 to make use of
%the FFT for wavelet transformation.
t_diff = diff(time); %find time interval
%The "secret sauce" equation below can be tuned by changing the divisor.
%Empirically 2 or 3 seems to be fine for pulses.
%However, increasing this value increases the wavelet sensitivity to noise.
%This is assuming the rate of p53 sampling is between 15 and 30 minutes.
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

function [out,wavelet_peaks] = findRidgeMap(in,s)
%Input:
%in: the CWT structure from the custom cwtft function in this file
%
%Output:
%out.scl: a cell array with ridges within each cell. Each ridge is a 1xN vector
%where N is the length the the ridge. The vector contains the scale component.
%out.cfs: a cell array the size of out.xy, but in each cell is a 1xN
%vector where N is the length of the ridge. The vector contains the scale
%of the ridge.
%out.t: a cell array the size of out.xy, but in each cell is a 1xN
%vector where N is the length of the ridge. The vector contains the time
%of the ridge.
%out.map: the ridge map, which is the same size of the wavelet transform
wavelet_peaks = cell(size(in.cfs,1),1);
%find smoothened signal
L = length(s);
s = s(1:round(L/10):end);
L2 = length(s);
s = smooth(s);
s = 0.01*interp1q((0:L2-1)',s,(0:(L2-1)/(L-1):L2-1)'); % 1% of the mean signal is the cutoff for meaningful peaks;
%find peaks at each scale
for i=1:size(in.cfs,1)
    wavelet_peaks{i} = first_pass_peak_detection(...
        in.cfs(i,:),in.scl(i)*2+1,s); %The window size was heursitically chosen to be (2*current_scale+1)
end
%-- Identify the Ridges --
ridge_map = zeros(size(in.cfs));
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
for i=size(in.cfs,1):-1:(gap_limit+1)
    for j=1:length(wavelet_peaks{i})
        for h=1:gap_limit
            %Search for peaks within the window size for scale i.
            low_bnd = wavelet_peaks{i}(j) - 2*in.scl(i-h);
            up_bnd = wavelet_peaks{i}(j) + 2*in.scl(i-h);
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
                    if (sum(wavelet_peaks_alias{i-h}==wavelet_peaks_alias{i}(j))==0) && (k_min<=k)
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
out.map = ridge_map; %just for show
[out.scl,out.cfs,out.t]=deal(cell(1,ridge_counter));
for i=1:ridge_counter
    for j=1:length(wavelet_peaks_alias)
        temp = (wavelet_peaks_alias{j} == i);
        if any(temp)
            if isempty(out.scl{i})
                out.scl{i}(1) = in.scl(j);
                out.sclind{i}(1) = j;
                out.t{i}(1) = wavelet_peaks{j}(temp);
                out.cfs{i}(1) = in.cfs(j,out.t{i}(1));
            else
                out.scl{i}(end+1)=in.scl(j);
                out.sclind{i}(end+1) = j;
                out.t{i}(end+1) = wavelet_peaks{j}(temp);
                out.cfs{i}(end+1) = in.cfs(j,out.t{i}(end));
            end
        end
    end
end
end

function [out]=first_pass_peak_detection(x,winw_size,s)
%There are many ways to find peaks. The wavelet method is powerful at
%detecting peaks but is actually dependent on simpler peak detection
%methods.

%Simple Peak Finding Method: Scan a signal with a window. Find the max. If the max is greater
%than the left and right endpoints of the window it is a peak candidate.
%The window is centered at each point of a waveform, so a wider peak
%candidate will recieve more votes in a sense. If a window has more than
%one point with the max value then the left most index is used. A downside
%to this method is that it is sensitive to the size of the window. However,
%this seeming limitation is actually put to use in the wavelet method by
%scaling the window along with the wavelet.
length_x = length(x);
if size(x,2) == 1
    x=x';
end
if winw_size<3
    winw_size = 3;
    warning('win_size:small','The window size for first pass peak detection may be too small and therefore yield illogical results.');
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
%Weed out unqualified peaks and peaks of questionable nature.
%If a peak has a negative wavelet coefficient value it is disqualified and
%if a peak is the first or last point of data it is disqualified.
for i=length(peak_elected):-1:1
    if x(peak_elected(i))<0
        peak_elected(i) = [];
    elseif peak_elected(i) == 1
        peak_elected(i) = [];
    elseif peak_elected(i) == length_x
        peak_elected(i) = [];
    end
end
%If a peak is less than the value of the 1% smoothened signal destroy its
%very existence. This is mainly to eliminate peaks that are very near to
%zero in the first few scales.
for i=length(peak_elected):-1:1
    if x(peak_elected(i))<s(peak_elected(i))
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

function []=rapid_report(wavelet_peaks,wavelet_xfrm_coefs,ridge_map,x,y)
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
penguinjet;
title('Ridge Map')
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
    pow102=2;
end

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
out.phz = 1./(out.scl*4); %0.25 is the pseudofrequency of the mexh for scale 1

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
    pow102=2;
end


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
out.phz = 1./(out.scl*5); %0.2 is the pseudofrequency of the dog1 for scale 1

out.cfs = real(out.cfs);
for i = 1:length(out.scl)
    out.cfs(i,:) = out.cfs(i,:)/sqrt(out.scl(i));
end
end

function [out,scale_stats] = processRidgeMap(in,scales,pseudoHz,len_sig)
%Input:
%in: the output, ridges, from findRidgeMap();
%
%Output:
%out: a struct with the peaks of the signal as determined by ridge analysis
%out.time = the time of a peak
%out.scale = the scale of that same peak
%out.waveletcfs = the wavelet coefficient at at that time and scale.
%scale_stats: statistics about the peaks contained in each scale. These are
%used to separate high frequency noise, from signal, from low frequency
%noise.
%scale_stats.numberOfPeaks = self explanatory
%scale_stats.meanCfsOfPeaks = mean value of the wavelet coefficients of the
%peaks at a given scale.
%scale_stats.varCfsOfPeaks = variance of the wavelet coefficients of the
%peaks at a given scale.
%scale_stats.peakEnrichment = Higher scales are expected to yield fewer
%peaks relative to lower scales. To account for this a measure of peak
%enrichment is defined as the number of peaks found relative to the number
%of peaks that would be identified if the signal was a sinusoid at the
%frequency that corresponds with that scale.

%pre-allocate struct that will contain peak information (assuming not more
%than 1000 peaks)
out(1000).time = [];
out(1000).scale = [];
out(1000).scaleindex = [];
out(1000).waveletcfs = [];
out(1000).ridgelength = [];
%Pre-allocate scale_stats
L = length(scales);
indArray = (1:L);
temp = num2cell(scales);
scale_stats = cell2struct(temp,'scale',1);
scale_stats(L).peakIdentity = [];
scale_stats(L).numberOfPeaks = [];
scale_stats(L).meanCfsOfPeaks = [];
scale_stats(L).varCfsOfPeaks = [];
scale_stats(L).peakEnrichment = [];
scale_stats(L).meanRidgeLength = [];
%Find the peak(s) of every ridge
ind = 1;
for i=1:length(in.cfs)
    %These next two lines of code represent an easy to implement solution
    %for finding peaks in 1D signals. This works especially well when the
    %signal is very smooth. The ridges should be very smooth given enough
    %scale density/coverage.
    peak_index = watershed(in.cfs{i});
    peak_index = find(~peak_index);
    temp_max = -inf;
    if ~isempty(peak_index)
        L = length(peak_index);
        for j=1:L
            out(ind).time = in.t{i}(peak_index(j));
            out(ind).scale = in.scl{i}(peak_index(j));
            out(ind).scaleindex = in.sclind{i}(peak_index(j));
            scale_stats(indArray(scales == out(ind).scale)).peakIdentity(end+1) = ind;
            out(ind).waveletcfs = in.cfs{i}(peak_index(j));
            
            if L>1
                switch(j)
                    case 1
                        out(ind).ridgelength = peak_index(j) + round(abs(peak_index(j+1)-peak_index(j))/2);
                    case L
                        out(ind).ridgelength = (length(in.cfs{i})-peak_index(j)) + round(abs(peak_index(j)-peak_index(j-1))/2);
                    otherwise
                        out(ind).ridgelength = round(abs(peak_index(j+1)-peak_index(j))/2) + round(abs(peak_index(j)-peak_index(j-1))/2);
                end
            else
                out(ind).ridgelength = length(in.cfs{i});
            end
            if temp_max < in.cfs{i}(peak_index(j))
                temp_max = in.cfs{i}(peak_index(j));
            end
            ind = ind+1;
        end
        if temp_max < max(in.cfs{i})
            [out(ind).waveletcfs,ind2] = max(in.cfs{i});
            out(ind).ridgelength = length(in.cfs{i});
            out(ind).time = in.t{i}(ind2);
            out(ind).scale = in.scl{i}(ind2);
            out(ind).scaleindex = in.sclind{i}(ind2);
            scale_stats(indArray(scales == out(ind).scale)).peakIdentity(end+1) = ind;
            ind = ind+1;
        end
    else
        [out(ind).waveletcfs,ind2] = max(in.cfs{i});
        out(ind).ridgelength = length(in.cfs{i});
        out(ind).time = in.t{i}(ind2);
        out(ind).scale = in.scl{i}(ind2);
        out(ind).scaleindex = in.sclind{i}(ind2);
        scale_stats(indArray(scales == out(ind).scale)).peakIdentity(end+1) = ind;
        ind = ind+1;
    end
end
out(ind:end)=[]; %remove empty pre-allocated space from struct

%Populate the scale statistics
pseudoHz = len_sig*pseudoHz;
for i=1:length(scale_stats)
    if isempty(scale_stats(i).peakIdentity)
        scale_stats(i).numberOfPeaks = 0;
        scale_stats(i).meanCfsOfPeaks = 0;
        scale_stats(i).varCfsOfPeaks = 0;
        scale_stats(i).peakEnrichment = 0;
        scale_stats(i).meanRidgeLength = 0;
    else
        scale_stats(i).numberOfPeaks = length(scale_stats(i).peakIdentity);
        temp = zeros(size(scale_stats(i).peakIdentity));
        temp2 = zeros(size(scale_stats(i).peakIdentity));
        for j=1:scale_stats(i).numberOfPeaks
            temp(j) = out(scale_stats(i).peakIdentity(j)).waveletcfs;
            temp2(j) = out(scale_stats(i).peakIdentity(j)).ridgelength;
        end
        scale_stats(i).meanCfsOfPeaks = mean(temp);
        scale_stats(i).varCfsOfPeaks = var(temp);
        scale_stats(i).peakEnrichment = scale_stats(i).numberOfPeaks/pseudoHz(i);
        scale_stats(i).meanRidgeLength = mean(temp2);
    end
end
end

function [out] = beautifyRidgeMap(pks,val,cfsmap)
%create the "white jet" colormap. It is the jet colormap, but the highest
%value
cfsmap_min = min(min(cfsmap));
cfsmap = ((cfsmap - cfsmap_min)*253/(max(max(cfsmap))-cfsmap_min))+1;
out = cfsmap;
out(pks>0) = 255;
out(val>0) = 0;
figure
penguinjet;
imagesc(out)
end

function [] = penguinjet()
tuxedojet = [0,0,0;0,0,0.53125;0,0,0.546875;0,0,0.5625;0,0,0.578125;0,0,0.59375;0,0,0.609375;0,0,0.625;0,0,0.640625;0,0,0.65625;0,0,0.671875;0,0,0.6875;0,0,0.703125;0,0,0.71875;0,0,0.734375;0,0,0.75;0,0,0.765625;0,0,0.78125;0,0,0.796875;0,0,0.8125;0,0,0.828125;0,0,0.84375;0,0,0.859375;0,0,0.875;0,0,0.890625;0,0,0.90625;0,0,0.921875;0,0,0.9375;0,0,0.953125;0,0,0.96875;0,0,0.984375;0,0,1;0,0.015625,1;0,0.03125,1;0,0.046875,1;0,0.0625,1;0,0.078125,1;0,0.09375,1;0,0.109375,1;0,0.125,1;0,0.140625,1;0,0.15625,1;0,0.171875,1;0,0.1875,1;0,0.203125,1;0,0.21875,1;0,0.234375,1;0,0.25,1;0,0.265625,1;0,0.28125,1;0,0.296875,1;0,0.3125,1;0,0.328125,1;0,0.34375,1;0,0.359375,1;0,0.375,1;0,0.390625,1;0,0.40625,1;0,0.421875,1;0,0.4375,1;0,0.453125,1;0,0.46875,1;0,0.484375,1;0,0.5,1;0,0.515625,1;0,0.53125,1;0,0.546875,1;0,0.5625,1;0,0.578125,1;0,0.59375,1;0,0.609375,1;0,0.625,1;0,0.640625,1;0,0.65625,1;0,0.671875,1;0,0.6875,1;0,0.703125,1;0,0.71875,1;0,0.734375,1;0,0.75,1;0,0.765625,1;0,0.78125,1;0,0.796875,1;0,0.8125,1;0,0.828125,1;0,0.84375,1;0,0.859375,1;0,0.875,1;0,0.890625,1;0,0.90625,1;0,0.921875,1;0,0.9375,1;0,0.953125,1;0,0.96875,1;0,0.984375,1;0,1,1;0.015625,1,0.984375;0.03125,1,0.96875;0.046875,1,0.953125;0.0625,1,0.9375;0.078125,1,0.921875;0.09375,1,0.90625;0.109375,1,0.890625;0.125,1,0.875;0.140625,1,0.859375;0.15625,1,0.84375;0.171875,1,0.828125;0.1875,1,0.8125;0.203125,1,0.796875;0.21875,1,0.78125;0.234375,1,0.765625;0.25,1,0.75;0.265625,1,0.734375;0.28125,1,0.71875;0.296875,1,0.703125;0.3125,1,0.6875;0.328125,1,0.671875;0.34375,1,0.65625;0.359375,1,0.640625;0.375,1,0.625;0.390625,1,0.609375;0.40625,1,0.59375;0.421875,1,0.578125;0.4375,1,0.5625;0.453125,1,0.546875;0.46875,1,0.53125;0.484375,1,0.515625;0.5,1,0.5;0.515625,1,0.484375;0.53125,1,0.46875;0.546875,1,0.453125;0.5625,1,0.4375;0.578125,1,0.421875;0.59375,1,0.40625;0.609375,1,0.390625;0.625,1,0.375;0.640625,1,0.359375;0.65625,1,0.34375;0.671875,1,0.328125;0.6875,1,0.3125;0.703125,1,0.296875;0.71875,1,0.28125;0.734375,1,0.265625;0.75,1,0.25;0.765625,1,0.234375;0.78125,1,0.21875;0.796875,1,0.203125;0.8125,1,0.1875;0.828125,1,0.171875;0.84375,1,0.15625;0.859375,1,0.140625;0.875,1,0.125;0.890625,1,0.109375;0.90625,1,0.09375;0.921875,1,0.078125;0.9375,1,0.0625;0.953125,1,0.046875;0.96875,1,0.03125;0.984375,1,0.015625;1,1,0;1,0.984375,0;1,0.96875,0;1,0.953125,0;1,0.9375,0;1,0.921875,0;1,0.90625,0;1,0.890625,0;1,0.875,0;1,0.859375,0;1,0.84375,0;1,0.828125,0;1,0.8125,0;1,0.796875,0;1,0.78125,0;1,0.765625,0;1,0.75,0;1,0.734375,0;1,0.71875,0;1,0.703125,0;1,0.6875,0;1,0.671875,0;1,0.65625,0;1,0.640625,0;1,0.625,0;1,0.609375,0;1,0.59375,0;1,0.578125,0;1,0.5625,0;1,0.546875,0;1,0.53125,0;1,0.515625,0;1,0.5,0;1,0.484375,0;1,0.46875,0;1,0.453125,0;1,0.4375,0;1,0.421875,0;1,0.40625,0;1,0.390625,0;1,0.375,0;1,0.359375,0;1,0.34375,0;1,0.328125,0;1,0.3125,0;1,0.296875,0;1,0.28125,0;1,0.265625,0;1,0.25,0;1,0.234375,0;1,0.21875,0;1,0.203125,0;1,0.1875,0;1,0.171875,0;1,0.15625,0;1,0.140625,0;1,0.125,0;1,0.109375,0;1,0.09375,0;1,0.078125,0;1,0.0625,0;1,0.046875,0;1,0.03125,0;1,0.015625,0;1,0,0;0.984375,0,0;0.96875,0,0;0.953125,0,0;0.9375,0,0;0.921875,0,0;0.90625,0,0;0.890625,0,0;0.875,0,0;0.859375,0,0;0.84375,0,0;0.828125,0,0;0.8125,0,0;0.796875,0,0;0.78125,0,0;0.765625,0,0;0.75,0,0;0.734375,0,0;0.71875,0,0;0.703125,0,0;0.6875,0,0;0.671875,0,0;0.65625,0,0;0.640625,0,0;0.625,0,0;0.609375,0,0;0.59375,0,0;0.578125,0,0;0.5625,0,0;0.546875,0,0;0.53125,0,0;0.515625,0,0;1,1,1];
colormap(tuxedojet);
end

function [ssas,f] = findSSAS(signal,time,sf)
%Input:
%signal: the original input signal
%time: the times at which the original signal was measured.
%sf: sampling frequency
%
%Output:
%ssas: the single-sided amplitude spectrum of the signal.

t_diff = diff(time); %find time interval
temp = median(t_diff);
t = (time(1):temp:time(end)); %time points are equally spaced
s = spline(time,signal,t); %Interpolate with splines
L = length(s);
h = hamming(L);
if size(h,1) == size(s,1)
    s = s.*h;
else
    s = s.*h';
end
NFFT = 2^nextpow2(L); % Next power of 2 from length of the signal, s
ssas = fft(s,NFFT)/L;
ssas = 2*abs(ssas(1:NFFT/2+1));
f = sf/2*linspace(0,1,NFFT/2+1);
% Plot single-sided amplitude spectrum.
% stem(f,ssas,'fill','--')
% title('Single-Sided Amplitude Spectrum')
% xlabel('Frequency (Hz)')
% ylabel('|Y(f)|')
end

function [outpks] = processPeaks(s,inpks,sts)
%initialize outpks
%Order the peak set first by scale and 2nd by time in ascending order. All
%peaks are given the same label of "type 0".
%type 1: high frequency noise
%type 2: signal
%type 3: low frequency noise
disp('here I am')

Lp=length(inpks);
Ls = length(s);
temp = (1:Lp)'; %Represents the original indexing found in inpks
ptime = [inpks(:).time];
pscl = [inpks(:).scale];
palltogether = [ptime',pscl',temp];
%First sort by scale
palltogether = sortrows(palltogether,2);
%Second sort by time
unique_scales = unique(pscl);
for i=1:length(unique_scales)
    temp1 = palltogether(:,2)==unique_scales(i);
    temp2 = palltogether(temp1,:);
    temp2 = sortrows(temp2,1);
    palltogether(temp1,:) = temp2;
end
%initialize the outpks
outpks = inpks(palltogether(:,3));
%set the type of all peaks to 0
outpks(end).type = [];
temp = zeros(1,Lp);
temp = num2cell(temp);
[outpks(:).type]=deal(temp{:});
%create a peak map
peakmap = zeros(length(sts),length(s));
for i=1:Lp
    j = outpks(i).scaleindex;
    k = outpks(i).time;
    peakmap(j,k) = i;
end
figure
imagesc(peakmap)
penguinjet;
%initialize the matchup matrix that keeps track of the pair of peaks where
%the current peak is less than the "greater" peak.
matchup = zeros(Lp,1);
%starting with the lowest and earliest peak search for nearby peaks in the
%area of a cone of influence
Lsc = length(sts);
stsScl = round(2.5*[sts(:).scale]); %The scale determines the size of the window of time
for i=1:Lp
    coi = zeros(Lsc,2);
    coi(:,1) = outpks(i).time - stsScl; %Left time point
    coi(:,2) = outpks(i).time + stsScl; %Right time point
    coi(coi<1) = 1;
    coi(coi>Ls) = Ls;
    stat0 = outpks(i).ridgelength*outpks(i).waveletcfs;
    flag1 = false; %flag1 ensures that the peak comparison stops once a "greater" peak is found.
    for j=outpks(i).scaleindex+1:Lsc
        if flag1
            break;
        end
        temp = peakmap(j,coi(j,1):coi(j,2));
        temp_cnt = sum(any(temp));
        switch(temp_cnt)
            case 0
                %do nothing
            case 1
                for k=1:length(temp)
                    if temp(k) ~= 0
                        stat1 = outpks(temp(k)).ridgelength*outpks(temp(k)).waveletcfs;
                        if stat0<stat1
                            if outpks(i).type == 0
                                outpks(i).type = 1;
                            end
                            matchup(i) = temp(k);
                            flag1 = true;
                        end
                        break;
                    end
                end
            otherwise
                coiDist = -inf;
                ind = 0;
                for k=1:length(temp)
                    if temp(k) ~= 0
                        tempDist = abs(outpks(i).time - outpks(temp(k)).time);
                        if tempDist > coiDist
                            ind = k;
                            coiDist = tempDist;
                        end
                    end
                end
                stat1 = outpks(temp(ind)).ridgelength*outpks(temp(ind)).waveletcfs;
                if stat0<stat1
                    if outpks(i).type == 0
                        outpks(i).type = 1;
                    end
                    matchup(i) = temp(k);
                    flag1 = true;
                end
        end
    end
    if ~flag1
        outpks(i).type = 2;
    end
end
%create a peak map
peakmap2 = zeros(length(sts),length(s));
for i=1:Lp
    j = outpks(i).scaleindex;
    k = outpks(i).time;
    peakmap2(j,k) = outpks(i).type;
end
figure
imagesc(peakmap2)
%Now scan all the type 2 peaks from the highest scale to the lowest. Look
%at the peaks in the lower scale in a block the width of the scale of the
%current peak. If there are any peaks below the current peak that is
%"greater" than the current peak becomes type 3.
%
%low frequency noise threshold: when two peaks occur very close together,
%their troughs will overlap. Think of two PSFs close to each other. To
%weight lower frequency peaks against this a threshold is used.
lfthresh = 1.25;
%isolate only type 2 peaks
temp = [outpks(:).type];
typeIndex = (1:Lp)';
temp = temp == 2;
typeIndex = typeIndex(temp)';
typeIndex = fliplr(typeIndex);
for i = typeIndex
    blk = zeros(1,2); %The block of time to search for peaks in
    stsScl = outpks(i).scale; %The scale determines the size of the window of time
    blk(1) = outpks(i).time - stsScl; %Left time point
    blk(2) = outpks(i).time + stsScl; %Right time point
    blk(blk<1) = 1;
    blk(blk>Ls) = Ls;
    stat0 = outpks(i).ridgelength*outpks(i).waveletcfs;
    flag1 = false; %flag1 ensures that the peak comparison stops once a "greater" peak is found.
    for j=outpks(i).scaleindex-1:-1:1
        if flag1
            break;
        end
        temp = peakmap(j,blk(1):blk(2));
        temp_cnt = sum(any(temp));
        switch(temp_cnt)
            case 0
                %do nothing
            case 1
                for k=1:length(temp)
                    if temp(k) ~= 0
                        stat1 = lfthresh*outpks(temp(k)).ridgelength*outpks(temp(k)).waveletcfs;
                        if stat0<stat1
                            outpks(i).type = 3;
                            flag1 = true;
                        end
                        break;
                    end
                end
            otherwise
                for k=1:length(temp)
                    if temp(k) ~= 0
                        stat1 = lfthresh*outpks(temp(k)).ridgelength*outpks(temp(k)).waveletcfs;
                        if stat0<stat1
                            outpks(i).type = 3;
                            flag1 = true;
                        end
                        break;
                    end
                end
        end
    end
end
%create a peak map
peakmap2 = zeros(length(sts),length(s));
for i=1:Lp
    j = outpks(i).scaleindex;
    k = outpks(i).time;
    peakmap2(j,k) = outpks(i).type;
end
figure
imagesc(peakmap2)
%Now revisit all the type 1 peaks and scan for the peaks in their cone of
%influence again. If the "greater" peak that made the current peak type 1
%in the first round is type 3, then make the current peak type 2.

%isolate only type 1 peaks
temp = [outpks(:).type];
typeIndex = (1:Lp)';
temp = temp == 1;
typeIndex = typeIndex(temp)';
for i=typeIndex
    if outpks(matchup(i)).type == 3
        outpks(i).type = 2;
    end
end
%create a peak map
peakmap2 = zeros(length(sts),length(s));
for i=1:Lp
    j = outpks(i).scaleindex;
    k = outpks(i).time;
    peakmap2(j,k) = outpks(i).type;
end
figure
imagesc(peakmap2)
%the signal peaks are adjusted to the highest nearby point on the signal.
%isolate only type 2 peaks
temp = [outpks(:).type];
typeIndex = (1:Lp)';
temp = temp == 2;
typeIndex = typeIndex(temp)';
typeIndex = fliplr(typeIndex);
for i = typeIndex
    time = outpks(i).time;
    window = outpks(i).scale;
    indL = time-window;
    if indL <= 0
        indL = 1;
    end
    indR = time+window;
    if indR > Ls
        indR = Ls;
    end
    [max_s,ind] = max(s(indL:indR));
    ind = ind + indL - 1;
    outpks(i).time = ind;
    outpks(i).value = max_s;
end
%plot peaks
figure
plot(s)
hold on
%Show signal peaks
temp1 = zeros(Ls,1);
temp2 = zeros(Ls,1);
for i=typeIndex
        temp1(i) = outpks(i).time;
        temp2(i) = outpks(i).value;
end
peak_index = temp1(temp1>0);
peak_value = temp2(temp2>0);
scatter(peak_index,peak_value,'filled','MarkerFaceColor','red');
end

function [] = plotPeaksAndValleys(s,pks,vly)
figure
plot(s)
hold on
%Show signal peaks
temp1 = zeros(length(s),1);
temp2 = zeros(length(s),1);
for i=1:length(pks)
    if pks(i).type == 2
        temp1(i) = pks(i).time;
        temp2(i) = pks(i).value;
    end
end
peak_index = temp1(temp1>0);
peak_value = temp2(temp2>0);
scatter(peak_index,peak_value,'filled','MarkerFaceColor','red');
%Show signal valleys
temp1 = zeros(length(s),1);
temp2 = zeros(length(s),1);
for i=1:length(vly)
    if vly(i).type == 2
        temp1(i) = vly(i).time;
        temp2(i) = vly(i).value;
    end
end
peak_index = temp1(temp1>0);
peak_value = temp2(temp2>0);
scatter(peak_index,peak_value,'filled','MarkerFaceColor','blue');
s2=s;
temp=(1:length(s));
temp(peak_index)=[];
s2(temp)=NaN;
plot(s2,'o','color','r')
hold off
end