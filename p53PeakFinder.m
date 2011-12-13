function [pks,vly,pksPitch,vlyPitch,sus,susPitch] = p53PeakFinder(signal,time,samplingFreq)


%Extend signal 25 percent both ends
%Interpolate the signal to increase the density of points 

%The wavelet transform is sensitive to edge effects. Since the signal
%(almost always) begins and ends abruptly, i.e. not equal to zero, the cone
%of influence from the edges will bleed into signal at high scales. To
%reduce this unwanted edge effect a tukey window is applied to the input
%signal. In addition, the first and last value of the signal are extended
%25 percent total signal length each direction for padding
%window the signal with a tukey window
tukey50=tukeywin(L,0.5);

%waveinfo('mexh')
%To locate both peaks, valleys, and their pitch use mexican hat wavelet
%waveinfo('haar')
%To locate singularities and switches use the haar wavelet.
end