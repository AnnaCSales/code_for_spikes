function [probe_fig] = plotProbe(spikeStruct)

% Plots positions of clusters on the probe.
probe_fig=figure('Name', 'Centre channel for each cluster', 'units', 'centimeters', 'pos', [2, 2, 12, 16]);
title('Position of cluster centre channels')
plot(spikeStruct.xcoords, spikeStruct.ycoords, 'or')
hold on
xlim([min(spikeStruct.xcoords)-80,  max(spikeStruct.xcoords)+80]);
n_cols=length(unique(spikeStruct.xcoords)); %the number of columns of recording sites on the probe.
ylim([min(spikeStruct.ycoords)-150,  max(spikeStruct.ycoords)+150]);

for chan=1:length(spikeStruct.xcoords)  %go through each channel, make a note of which clusters are there.
    clusts_here=find(spikeStruct.c_channel==chan);
    if clusts_here
        if mod(chan, n_cols)==0
            text(spikeStruct.xcoords(chan)+40,spikeStruct.ycoords(chan),[num2str(clusts_here)],'HorizontalAlignment','right')
        elseif mod(chan,n_cols)==n_cols-1;
            text(spikeStruct.xcoords(chan)-40,spikeStruct.ycoords(chan),[ num2str(clusts_here)],'HorizontalAlignment','left')
        end
    end
end

ax1 = gca;                   % gca = get current axis
ax1.YAxis.Visible = 'off';   % remove y-axis
ax1.XAxis.Visible = 'off';   % remove y-axis
%draw a probe outline
xs=[min(spikeStruct.xcoords)-10, min(spikeStruct.xcoords)-10, mean(spikeStruct.xcoords),max(spikeStruct.xcoords)+10,max(spikeStruct.xcoords)+10];
ys=[max(spikeStruct.ycoords)+20, min(spikeStruct.ycoords), min(spikeStruct.ycoords)-80,min(spikeStruct.ycoords), max(spikeStruct.ycoords)+20];
patch(xs, ys, 'k', 'FaceColor','blue','FaceAlpha',.2)
title('Position of cluster centre channels')

end

