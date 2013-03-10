function [myConfig] = openConfig()
%Find the root path
myConfig = my_parseXML('config.m');
myConfig = myConfig.p53Cinema_config;
%Check for existence of output directory
if ~exist(myConfig.output_directory.txtd8a{1},'dir')
    error('p53CinemaConfig:invalidOutputDir','An invalid output directory was found\nwhen opening the configuration file.');
end