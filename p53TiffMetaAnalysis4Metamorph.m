function [time] = p53TiffMetaAnalysis4Metamorph(t)
% time = p53TiffMetaAnalysis4Metamorph(t)
% Input:
% t = a TIFF image object created by "Tiff(filename,'r')".
% Output:
% time = the time the image was acquired by the microscope. 'YYYYMMDD
% 24:60:60.000'
% Description:
% The TIFF format can contain useful metadata about the image. When a stack
% of images is created the metadata from the individual images will be
% consolidated into the "ImageDescription" TIFF tag in XML format within the
% first directory of the TIFF stack. Metamorph stores a copious amount of
% info about the experimental setup and the time of acquisition. 
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
% -----  -----
metadata = t.getTag('ImageDescription');
fid = fopen('t3mp.xml','w');
fprintf(fid,'%s',metadata);
fclose(fid);
xdoc = xmlread('t3mp.xml');
node_root = xdoc.getDocumentElement;
my_list = node_root.getElementsByTagName('prop');
for i=0:my_list.getLength-1
    if my_list.item(i).hasAttributes
        if strcmp(my_list.item(i).getAttribute('id').toString.toCharArray','acquisition-time-local')
            time = my_list.item(i).getAttribute('value').toString.toCharArray';
        end
    end
end
delete('t3mp.xml');
end