function []=flatfieldcorrection(stackpath,ffpath)
% [] = flatfieldcorrection(stackpath,ffpath)
% Input:
% stackpath: the location of TIFF that will be flat field corrected
% ffpath: the location of the correction images 
%
% Output:
% There is no direct argument output. Rather, new stacks will be created from
% the stacks found in the stackpath directory and stored in a folder called
% [stackpath '_ff'].
%
% Description:
% 
%
% Other Notes:
% 
if nargin==0
    stackpath = 'C:\Users\Kyle\Documents\p53CinemaSTACKS';
    ffpath='G:\KWKDocuments\My Dropbox\p53Cinema\flatfield_20110909';
end
%----- Import stacknames -----
dirCon_stack = dir(stackpath);
disp(['working in ', stackpath]); %Sanity Check
expr='.*_w\d+.*_s\d+.*_t.*\.tif';
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
stacknames=Temp;
dirCon_ff = dir(ffpath);
end

