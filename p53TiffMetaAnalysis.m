function [] = p53TiffMetaAnalysis(t)
% S = my_parseXML(filename)
% Input:
% t = a TIFF image object created by "Tiff(filename,'r')".
% Output:
%
% Description:
% The TIFF format can contain useful metadata about the image. When a stack
% of images is created the metadata from the individual images will be
% consolidated into the "ImageDescription" TIFF tag in XML format
% Other Notes:
% 
                %<parsing the tiff file with these commands>
                %             imfinfo(name2read,'tif');
                %             t = Tiff(filename,'r');
                %             t.nextDirectory;
                %             t.setDirectory;
                %             t.lastDirectory;
                %             t.currentDirectory;
                %             t.writeDirectory;
                %             t.rewriteDirectory;
                %             t.getTag('ImageLength');
                %<
                
end