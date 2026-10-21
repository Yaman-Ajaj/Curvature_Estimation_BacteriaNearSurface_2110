%% reset
clear all, close all, clc;

%% load data
load('Ref0722_unfiltered.mat')
load('Ref0722_filtered.mat')

%% 3d plot
figure(9901)
for ii = 1 : 300
    plot3(AllBugs_unfilt{ii}(:,2),AllBugs_unfilt{ii}(:,3),AllBugs_unfilt{ii}(:,4),'-o')
    hold on
    plot3(AllBugs_filt{ii}(:,2),AllBugs_filt{ii}(:,3),AllBugs_filt{ii}(:,4),'-o')
    hold off
    title(num2str(ii))
    
    % press 'n' to continue
    h = annotation('textbox', [0.65, 0.02, 0.3, 0.05],'String','Press "n" to continue','EdgeColor', 'none','HorizontalAlignment','right');
    set(gcf, 'KeyPressFcn', @(src, event) setappdata(gcf, 'key', event.Key));
    setappdata(gcf, 'key', '');
        while true
            key = getappdata(gcf, 'key');
            if strcmp(key, 'n')
                break
            end
            pause(0.05)
        end
    clf;
end

