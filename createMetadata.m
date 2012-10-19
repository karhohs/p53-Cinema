docNode = com.mathworks.xml.XMLUtils.createDocument('metadata');
metadata = docNode.getDocumentElement;
%time of Acquisition
timeOfAcquisition = docNode.createElement('timeOfAcquisition');
timeOfAcquisitionText = '1293';
timeOfAcquisition.appendChild(docNode.createTextNode(timeOfAcquisitionText));
metadata.appendChild(timeOfAcquisition);
%timepoint
timepoint = docNode.createElement('timepoint');
timepointText = '3';
timepoint.appendChild(docNode.createTextNode(timepointText));
metadata.appendChild(timepoint);
%wavelength
wavelength = docNode.createElement('wavelength');
wavelengthType = docNode.createElement('type');
wavelengthTypeText = 'YFP';
wavelengthType.appendChild(docNode.createTextNode(wavelengthTypeText));
wavelength.appendChild(wavelengthType);
wavelengthExposure = docNode.createElement('exposure');
wavelengthExposureText = '100ms';
wavelengthExposure.appendChild(docNode.createTextNode(wavelengthExposureText));
wavelength.appendChild(wavelengthExposure);
metadata.appendChild(wavelength);
%position
position = docNode.createElement('position');
positionText = '2';
position.appendChild(docNode.createTextNode(positionText));
metadata.appendChild(position);
%z-slice
zSlice = docNode.createElement('z-slice');
zSliceText = '1';
zSlice.appendChild(docNode.createTextNode(zSliceText));
metadata.appendChild(zSlice);
%label
label = docNode.createElement('label');
labelText = 'mcf7';
label.appendChild(docNode.createTextNode(labelText));
metadata.appendChild(label);
%write to a file
filename = 'test.xml';
xmlwrite(filename,docNode);