function [] = CloneSelectCP(path2mat)
%Find all the file names of the images analyzed by Cell Profiler.
%Use the corresponding CloneSelectCP.cp pipeline to generate the data.
data = open(path2mat);
%Find all the unique well names and how many images were taken in each
%well.
filenames = data.handles.Measurements.Image.FileName_p53;
wellNamesAll = cell(size(filenames));
for i = 1:length(filenames)
    wellNamesAll{i} = regexp(filenames{i},'(?<=_).*(?=_s)','match','once');
end
[wellNamesUnique, ~, ind] = unique(wellNamesAll);
numberOfImagesPerWell = histc(ind, 1:length(wellNamesUnique));
%Find the indices of all the images taken in each well from the filename
%cell
filenameIndex = cell(1,length(wellNamesUnique));
filenameIndicies = (1:length(filenames));
for i = 1:length(wellNamesUnique)
    my_counter = 0;
    indexEntry = zeros(1,numberOfImagesPerWell(i));
    for j = filenameIndicies
        my_match = regexp(filenames{j},wellNamesUnique{i}, 'once');
        if ~isempty(my_match)
            my_counter = my_counter + 1;
            indexEntry(my_counter) = j;
        end
        if my_counter == numberOfImagesPerWell(i)
            break
        end
    end
    filenameIndex{i} = indexEntry;
end
%For each well and all the images taken in that well...
cloneSelectStats(length(wellNamesUnique)) = struct('mean',[],'cv',[],'count',[],'well',[]);
for i = 1:length(wellNamesUnique)
    indexEntry = filenameIndex{i};
    %Find the mean intensity for the cell population and the number of
    %cells identified and the coefficient of variation for the cell
    %population.
    temp = data.handles.Measurements.Nuclei.('Intensity_IntegratedIntensity_p53 corrected');
    integratedIntensities = [];
    for j = indexEntry
        integratedIntensities = cat(1,integratedIntensities,temp{j});
    end
    cloneSelectStats(i).mean = mean(integratedIntensities);
    cloneSelectStats(i).cv = std(integratedIntensities)/cloneSelectStats(i).mean;
    cloneSelectStats(i).count = numel(integratedIntensities);
    %Identify the name of the well
    cloneSelectStats(i).well = wellNamesUnique{i};
end

%Calculate a few rank statistics (We assume most wells are not empty)
cloneMean = [cloneSelectStats(:).mean];
realCloneCV = [cloneSelectStats(:).cv];
cloneCV = 1./realCloneCV;
cloneCount = [cloneSelectStats(:).count];
%Screen for empty wells using mean intensity to separate
my_iqr = iqr(cloneMean);
my_std = 0.7413*my_iqr; %har
my_mean = mean(cloneMean);
my_cutoff = my_mean - 2.5*my_std;
goodIndex = cloneMean>my_cutoff;
cloneMean = cloneMean(goodIndex);
cloneCV = cloneCV(goodIndex);
cloneCount = cloneCount(goodIndex);
wellNamesUnique = wellNamesUnique(goodIndex);
%Rank each well by the mean, cv, and number of cells
ind = 1:length(wellNamesUnique);
rankMean = [ind;cloneMean]';
rankMean = sortrows(rankMean,2);
rankInd = [ind',rankMean(:,1)];
rankInd = sortrows(rankInd,2);
rankMean = rankInd(:,1);
rankCV = [ind;cloneCV]';
rankCV = sortrows(rankCV,2);
rankInd = [ind',rankCV(:,1)];
rankInd = sortrows(rankInd,2);
rankCV = rankInd(:,1);
rankCount = [ind;cloneCount]';
rankCount = sortrows(rankCount,2);
rankInd = [ind',rankCount(:,1)];
rankInd = sortrows(rankInd,2);
rankCount = rankInd(:,1);
%Create a cumulative rank by adding the ranks from mean, cv, and number of
%cells together. The ranks can be weighted
rank = rankMean + 0.5*rankCV + 0.5*rankCount;
rankInd = [ind',rank];
rank = sortrows(rankInd,2)';
%Create a plot that summarizes the findings
S = cloneCV(rank(1,:))/max(cloneCV(rank(1,:)))*20;
S = S.^2;
colormap(hot)
scatter(ind,cloneMean(rank(1,:)),S,cloneCount(rank(1,:)),'filled')
set(gca,'FontName','Arial','FontSize',12,'XTick',(1:round(ind(end)/12):ind(end)),'XTickLabel',wellNamesUnique(rank(1,end-12:end)))
end