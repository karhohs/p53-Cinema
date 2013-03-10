function [] = p53plot(path2mat,folderout)
global fout
fout = folderout;
data = open(path2mat);
signal = cell(length(data.unitOfLife),1);
[signal{:}] = deal(data.unitOfLife(:).meanIntensity);
time = cell(length(data.unitOfLife),1);
[time{:}] = deal(data.unitOfLife(:).time);
traces(signal,time)
clear global fout
end

function [] = traces(levels,time)
global fout
h=figure ('visible', 'off', 'position', [10, 10, 1920, 1080]);
ax=axes('parent',h);
set(ax,'fontsize',30);
plot(time{1},levels{1})
%xlim(ax,[0 100]);
xlabel('time (hours)','fontsize',34);
%ylim(ax,[0 0.12]);
ylabel('p53 Signal (A.U.)','fontsize',34);
title('Mean Fluorescent Intensity of p53 Signal','fontsize',34);
%hold on
%hold off
%saveas (ax, filename, 'pdf' );
print(h,fullfile(fout,'traces.eps'),'-dpsc2','-painters');
close(h);
if length(time)>1
    for i=2:length(time)
        h=figure ('visible', 'off');
        ax=axes('parent',h);
        set(ax,'fontsize',30);
        plot(time{i},levels{i})
        %xlim(ax,[0 100]);
        xlabel('time (hours)');
        %ylim(ax,[0 0.12]);
        ylabel('p53 Signal (A.U.)');
        title('Mean Fluorescent Intensity of p53 Signal');
        %hold on
        %hold off
        %saveas (ax, filename, 'pdf' );
        print(h,fullfile(fout,'traces.eps'),'-dpsc2','-painters','-append');
        close(h);
    end
end
end