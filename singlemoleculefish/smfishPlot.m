function [] = smfishPlot(path2mat,folderout,path2maxproj,fn)
global fout
global maxproj
global filename
maxproj = imread(path2maxproj,'tiff');
fout = folderout;
filename = fn;
data = open(path2mat);
rings(data.fociarray,data.sizeOfImage)
clear global fout
clear global maxproj
clear global filename
end

function [] = rings(fociarray,sizeOfImage)
global fout
global maxproj
global filename
linewidth = 1;
circleradiusbase = 5; %units in pixels
circlecircumference = 2*pi*circleradiusbase;
ringIM = zeros(sizeOfImage(1:2));
for k = 1:linewidth;
    if linewidth > 1
        circleradius = circleradiusbase + k - linewidth/2;
    else
        circleradius = circleradiusbase;
    end
    for i = 1:length(fociarray)
        [y2,x2,~] = ind2sub(sizeOfImage,fociarray(i));
        for j=1:round(2*circlecircumference)
            circlex = round(x2 + circleradius*cos(j*2*pi/round(2*circlecircumference)));
            circley = round(y2 + circleradius*sin(j*2*pi/round(2*circlecircumference)));
            if circlex < 1
                continue
            elseif circlex > sizeOfImage(2)
                continue
            end
            if circley < 1
                continue
            elseif circley > sizeOfImage(1)
                continue
            end
            ringIM(circley,circlex) = 1;
        end
    end
end
ringIM = logical(ringIM);
ringPicture = zeros(sizeOfImage(1),sizeOfImage(2),3);
ringPicture(:,:,1) = 255*ringIM;
ringPicture(:,:,2) = 255*ringIM;
maxProj2 = zeros(sizeOfImage(1),sizeOfImage(2),3);
mp1 = maxproj;
mp1(ringIM) = 255;
mp2 = maxproj;
mp2(ringIM) = 255;
mp3 = maxproj;
mp3(ringIM) = 0;
maxProj2(:,:,1) = mp1;
maxProj2(:,:,2) = mp2;
maxProj2(:,:,3) = mp3;
name = regexprep(filename,'(\w*)(?=\.)','$1_circlesANDmaxproj');
imwrite(uint8(maxProj2),fullfile(fout,name),'tif','Compression','none');
name = regexprep(filename,'(\w*)(?=\.)','$1_circles');
imwrite(uint8(ringPicture),fullfile(fout,name),'tif','Compression','none');
end

function [] = graphtemplate()
h=figure ('visible', 'off', 'position', [10, 10, 1920, 1080]);
ax=axes('parent',h);
set(ax,'fontsize',30);
xlim(ax,[0 100]);
xlabel('Total Fluorescent Intensity of DNA Signal (A.U.)','fontsize',34);
ylim(ax,[0 0.12]);
ylabel('Probability Density','fontsize',34);
title('Total Fluorescent Intensity of p53 Signal (A.U.), 0hr','fontsize',34);
hold on
hold off
saveas (ax, filename, 'pdf' );
print(h,fullfile(fout,'new_image.eps'),'-depsc2');
close(h);
end