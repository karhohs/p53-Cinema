function [] = lowCopyNumberThresholdFinder(smfishpath,varargin)
    % [] = my_tiffStacker(path,positions,timepoints)
% Input:
% path: a char. The path to the folder that contains the raw TIFF images from
% imageJ.
% positions (optional): an array of positive integers. The integers in this array represent the
% positions for which stacks will be created.
% timepoints (optional): a cell containing arrays of positive integers or
% an array of positive integers. If it is a cell, then it must contain an
% array for each position. If it is an array, only these timepoints will be
% taken for each position.
%
% Output:
% There is no direct argument output. Rather, a new folder will be created
% that recomputes the content of the smfishStackPreprocessing
%
% Description:
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