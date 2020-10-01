function [wavefig] = plot_across_chans(spikeStruct, unit)
%Plots waveforms across multiple channels around the centre channel

cent_chan=spikeStruct.c_channel(unit);
startchan=cent_chan-2;
endchan=cent_chan+2;

y_plt=0.1;
left_val=0.2;
right_val=0.65;
if mod(startchan,2) %if it's odd, it's on LHS, lower
    x_plt=left_val;
else
    x_plt=right_val;
end

wfs=squeeze(spikeStruct.allchanWFs(unit, :, :));
wfs_sem=squeeze(spikeStruct.allchanSTDs(unit, :,:));
wfs_n_extracted=spikeStruct.nWFs_extracted(unit);

wavefig=figure('Color', 'w')
wave_time=1000/spikeStruct.sample_rate *(1:size(wfs,2)); %convert to milliseconds
   % plot on centre channel, two units above, two units below.
chans_to_plot=startchan:endchan;

if startchan<3
    chans_to_plot(chans_to_plot<1)=[];
    chans_to_plot=[chans_to_plot, chans_to_plot(end)+1,chans_to_plot(end)+2 ];
end


ybounds=round(max(abs(wfs(cent_chan,:))),-1)+10;
for pl=1:5
    chan_=chans_to_plot(pl);
    if chan_==cent_chan
        pltcol='r';
    else
        pltcol='b';
    end
   
    
    
    subplot('Position', [x_plt, y_plt, 0.2, 0.18])
    shadedErrorBarLight(wave_time, wfs(chan_, :), wfs_sem(chan_,:)./sqrt(wfs_n_extracted), pltcol, 1)
    xlabel('ms');
    ylabel('\mu V');
   
    ylim([-ybounds, ybounds]);
    xlim([0, wave_time(end)]);
    if x_plt==left_val
        x_plt=right_val;
    else
        x_plt=left_val;
    end
    y_plt=y_plt+.15;
    title(['Chan: ' num2str(chan_) ' (OEP chan: ' num2str(1+spikeStruct.chanMap(chan_)) '.)'], 'Fontweight', 'normal') 
    ax = gca;
    ax.FontSize = 8;
   
end

if x_plt==left_val
    textx_pos=-5;
else
    textx_pos=6;
end

text(textx_pos, ybounds+30, ['Unit: ' num2str(unit)]);