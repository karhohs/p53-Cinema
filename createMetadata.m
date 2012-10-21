function [] = createMetadata(p)
docNode = com.mathworks.xml.XMLUtils.createDocument('metadata');
metadata = docNode.getDocumentElement;
%time of Acquisition
timeOfAcquisition = docNode.createElement('timeOfAcquisition');
timeOfAcquisition.appendChild(docNode.createTextNode(p.timeOfAcquisitionText));
metadata.appendChild(timeOfAcquisition);
%timepoint
timepoint = docNode.createElement('timepoint');
timepoint.appendChild(docNode.createTextNode(p.timepointText));
metadata.appendChild(timepoint);
%wavelength
wavelength = docNode.createElement('wavelength');
wavelengthType = docNode.createElement('type');
wavelengthType.appendChild(docNode.createTextNode(p.wavelengthTypeText));
wavelength.appendChild(wavelengthType);
% wavelengthExposure = docNode.createElement('exposure');
% wavelengthExposure.appendChild(docNode.createTextNode(p.wavelengthExposureText));
% wavelength.appendChild(wavelengthExposure);
metadata.appendChild(wavelength);
%position
position = docNode.createElement('position');
position.appendChild(docNode.createTextNode(p.positionText));
metadata.appendChild(position);
%z-slice
zSlice = docNode.createElement('z-slice');
zSlice.appendChild(docNode.createTextNode(p.zSliceText));
metadata.appendChild(zSlice);
%label
label = docNode.createElement('label');
label.appendChild(docNode.createTextNode(p.labelText));
metadata.appendChild(label);
%write to a file
xmlwrite(p.filename,docNode);