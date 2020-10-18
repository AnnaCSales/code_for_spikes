function [wavefig] = plot_across_chans(spikeStruct, unit, dist)
%Plots waveforms across multiple channels around the centre channel
%dist is the range in um for channels, from the centre channel
%Get waveforms out of spikeStruct
wfs=squeeze(spikeStruct.allchanWFs(unit, :, :));
wfs_sem=squeeze(spikeStruct.allchanSTDs(unit, :,:));
wfs_n_extracted=spikeStruct.nWFs_extracted(unit);
wave_time=1000/spikeStruct.sample_rate *(1:size(wfs,2)); %convert to milliseconds

%now work out where to plot, to accruately represent channel layout
cent_chan=spikeStruct.c_channel(unit);  %the centre channel for this unit

xs=spikeStruct.xcoords;
ys=spikeStruct.ycoords;


%find all channels within 100um of the centre
cent_y=ys(cent_chan);
chans_within_range=find(ys>cent_y-dist & ys<cent_y+dist);

%now scale channel coords from 0.1 to 0.75, to provide coordinates for the axes
%representing each channel
plotPos=[rescale(xs(chans_within_range),0.15, 0.75), rescale(ys(chans_within_range), 0.06, 0.85)];

%work out which channel is in the bottom left of the plot - this is where
%we will show axis info

[low_val,low_ind]=min(plotPos(:,2));
lowest=find(plotPos(:,2)==low_val);
[~,leftest]=min(plotPos(lowest, 1));
bottom_left_chan=chans_within_range(lowest(leftest));


wavefig=figure('Color', 'w', 'Units', 'normalized', 'Position', [0.1 0.1 0.4 0.75]);
for sb=1:length(chans_within_range)
    chan_=chans_within_range(sb);
    
    if chan_==cent_chan  %plot the centre channel in red to highlight
        pltcol='r';
    else
        pltcol='b';
    end
   
    %put the plot in the correct place.
    chanPlots(sb)=subplot('Position', [plotPos(sb,1), plotPos(sb,2), 0.18, 0.1]);
    shadedErrorBarLight(wave_time, wfs(chan_, :), wfs_sem(chan_,:)./sqrt(wfs_n_extracted), pltcol, 1);
    ylim([-150, 150]);
    xlim([0,wave_time(end) ]);
    title(['Chan: ' num2str(chan_) ' (OEP chan: ' num2str(1+spikeStruct.chanMap(chan_)) '.)'], 'Fontweight', 'normal') ;
    ax = gca;
    ax.FontSize = 8;
    
    
    if chan_==bottom_left_chan
        xlabel('ms');
        ylabel('\mu V');
        box off
    else
       set(gca,'xtick',[])
       set(gca,'ytick',[])
       box off
       axis off
    end
 
    
end
subplot(chanPlots(1))
aa=gca;
text(3.3, 150, ['Unit: ' num2str(unit)]);
text(3.3, 0, ['Centre channel in red.']);
text(3.3, -80, ['Showing channels within ' num2str(dist) '\mum']);