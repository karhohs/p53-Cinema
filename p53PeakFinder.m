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
% 

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

%Baseline removal: Remove the baseline using the discrete wavelet transform

%Find the ridgemap for peak detection using the continuous wavelet
%transform.


end