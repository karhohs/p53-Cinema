function [] = exportCellData2XML(uniqueUOL,outpath)
docNode = com.mathworks.xml.XMLUtils.createDocument('data');
data = docNode.getDocumentElement;
popStat = docNode.createElement('population_statistics');
%create an array that contains every signal measurement
allSignal = [];
for i=1:length(uniqueUOL)
    allSignal = cat(2,allSignal,uniqueUOL(i).meanIntensity');
end
popStat.setAttribute('max',num2str(max(allSignal)));
popStat.setAttribute('min',num2str(min(allSignal)));
popStat.setAttribute('mean',num2str(mean(allSignal)));
popStat.setAttribute('std',num2str(std(allSignal)));
data.appendChild(popStat);
%fill the XML file with signal data for each cell
for i=1:length(uniqueUOL)
    newCell = docNode.createElement('cell');
    newCell.setAttribute('id',num2str(i));
    data.appendChild(newCell);
    %add signal data for mean_intensity
    signal = docNode.createElement('signal');
    signal.setAttribute('type','mean_intensity');
    %make data comma delimited
    str = sprintf('%0.5e,',uniqueUOL(i).meanIntensity);
    str(end) = []; %remove the extra comma at the end
    signal.appendChild(docNode.createTextNode(str));
    newCell.appendChild(signal);
end
xmlwrite(fullfile(outpath,'data.xml'),docNode);
end