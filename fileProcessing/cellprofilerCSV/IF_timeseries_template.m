function [] = IF_timeseries_template(path)
if isempty(path)
    path = 'M:\Working\SI_201211B_ff\A549\48h';
end
startpath = pwd;
%% ----- Import PNG paths -----
%Check whether the computer is a mac
if ismac
    error('manSeg:notApc','Flatfield correction is currently configured to run on a pc');
    %[~,filepaths]=system('ls -Rp */*.png');
elseif ispc
    %the 'dir' command in Matlab does not search subfolders. However,
    %there is a MS-DOS command that will do the trick!
    %'/B' uses bare format (no heading information or summary)
    %'/S' displays files in specified directory and all subdirectories
    cd(path)
    [~,filepaths] = system('dir /S /B *.png');
    filepaths = textscan(filepaths,'%s','Whitespace','\b\t');
    filepaths = filepaths{1};
else
    error('manSeg:notApc','IF_timeseries_template.m is currently configured to run on a pc');
end
header = {'Image_FileName_Cy5',... %1
    'Image_PathName_Cy5',... %2
    'Image_FileName_FITC',... %3
    'Image_PathName_FITC',... %4
    'Image_FileName_DAPI',... %5
    'Image_PathName_DAPI'}; %6
M = cell(10000,length(header));
for i=1:length(filepaths)
    [pathstr,fname,ext] = fileparts(filepaths{i});
    reg = regexp(fname,'(?<=_w\d+)[^_]+','match');
    if isempty(reg)
        continue;
    end
    reg = reg{1};
    regsite = regexp(fname,'(?<=_s)\d+','match');
    regsite = str2double(regsite{1});
    switch reg
        case 'Cy5'
            M{regsite,1} = strcat(fname,ext);
            M{regsite,2} = pathstr;
        case 'FITC'
            M{regsite,3} = strcat(fname,ext);
            M{regsite,4} = pathstr;
        case 'DAPI'
            M{regsite,5} = strcat(fname,ext);
            M{regsite,6} = pathstr;
    end
end
isemptyM = cellfun(@isempty,M(:,1));
M(isemptyM,:) = [];
M = vertcat(header,M);
cd(startpath);
fid=fopen('cpCSV.csv','w');
csvFun = @(str)sprintf('%s,',str);
for i=1:size(M,1)
xchar = cellfun(csvFun, M(i,:), 'UniformOutput', false);
xchar = strcat(xchar{:});
fprintf(fid,'%s\n',xchar(1:end-1));
end
fclose(fid);
end