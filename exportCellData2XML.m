function [] = exportCellData2XML(uniqueUOL,outpath)
docNode = com.mathworks.xml.XMLUtils.createDocument('data');
data = docNode.getDocumentElement;
popStat = docNode.createElement('population_statistics');
%create an array that contains every signal measurement
allSignal = [];
allTime = [];
for i=1:length(uniqueUOL)
    allSignal = cat(2,allSignal,uniqueUOL(i).meanIntensity');
    allTime = cat(2,allTime,uniqueUOL(i).timePoints);
end
popStat.setAttribute('maxSignal',num2str(max(allSignal)));
popStat.setAttribute('minSignal',num2str(min(allSignal)));
popStat.setAttribute('meanSignal',num2str(mean(allSignal)));
popStat.setAttribute('stdSignal',num2str(std(allSignal)));
popStat.setAttribute('maxTime',num2str(max(allTime)));
popStat.setAttribute('minTime',num2str(max(allTime)));
data.appendChild(popStat);
%fill the XML file with signal data for each cell
for i=1:length(uniqueUOL)
    newCell = docNode.createElement('cell');
    newCell.setAttribute('id',num2str(i));
    data.appendChild(newCell);
    %-----add signal data for mean_intensity
    signal = docNode.createElement('signal');
    signal.setAttribute('type','mean_intensity');
    %make data comma delimited
    str = sprintf('%0.5e,',uniqueUOL(i).meanIntensity);
    str(end) = []; %remove the extra comma at the end
    signal.appendChild(docNode.createTextNode(str));
    newCell.appendChild(signal);
    %-----add data about divisions
    divisionTimes = docNode.createElement('division_times');
    divisionTimes.setAttribute('number_of_divisions',num2str(uniqueUOL(i).divisionbool));
    %make data comma delimited
    str = sprintf('%d,',uniqueUOL(i).divisionTime);
    str(end) = []; %remove the extra comma at the end
    divisionTimes.appendChild(docNode.createTextNode(str));
    newCell.appendChild(divisionTimes);
    %-----add data about death
    deathTime = docNode.createElement('death_time');
    deathTime.setAttribute('deathbool',num2str(uniqueUOL(i).deathbool));
    %make data comma delimited
    str = sprintf('%d,',uniqueUOL(i).deathTime);
    str(end) = []; %remove the extra comma at the end
    deathTime.appendChild(docNode.createTextNode(str));
    newCell.appendChild(deathTime);
    %-----add data about timepoints
    timePoints = docNode.createElement('time_points');
    %make data comma delimited
    str = sprintf('%d,',uniqueUOL(i).timePoints);
    str(end) = []; %remove the extra comma at the end
    timePoints.appendChild(docNode.createTextNode(str));
    newCell.appendChild(timePoints);
end
xmlwrite(fullfile(outpath,'data.xml'),docNode);
end