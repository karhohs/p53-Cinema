function [] = lowCopyNumberThresholdFinder(smfishpath,varargin)
    % [] = my_tiffStacker(path,positions,timepoints)
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
%----- flatfield correction for each z-slice -----
if p.Results.manualthresh

end

end