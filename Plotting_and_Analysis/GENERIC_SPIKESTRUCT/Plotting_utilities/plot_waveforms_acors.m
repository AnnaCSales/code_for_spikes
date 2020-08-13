function [wave_fig, auto_fig] = plot_waveforms_acors(spikeStruct, params)
% A plotting script that takes with a spikestruct, plots:

%  - waveforms
%  - autocors

% Dependencies: gausskernel, shadedErrorBar

% Expects to be passsed (a) the spikeStruct (b) a params stuct, with the
% following fields

% params.cell_list  - a vector containing the cluster numbers to plot.
% params.binsize - binsize, in ms
% params.window - window over which to plot the auto-cor.
% params.baseline - 1 if we're plotting only over baseline, 0 if not.

% Updated July 2020, Anna Sales.

%% Specify the subset of cells to plot, and the time period to consider.

GDcells=params.cell_list;

if params.baseline==1
    min_t=spikeStruct.baseline_st;
    max_t=spikeStruct.baseline_end;
else
    min_t=0;
    max_t=spikeStruct.timeRange(2);
end

xcor_window=params.window/1000;
%% Some other useful info from the spikeStruct.

nclusts=length(GDcells);
fs=spikeStruct.sample_rate;
newWFs=spikeStruct.allchanWFs;
av_waveform=spikeStruct.av_waveform;  %av waveform on the centre channel
bl_start=spikeStruct.baseline_st;  %define baseline, if needed
bl_end=spikeStruct.baseline_end;
min_t=spikeStruct.timeRange(1);  %define period of entire recording
max_t=spikeStruct.timeRange(2);
sq_=ceil(length(GDcells)^0.5) ; %for calculating the number of subplots required

%% Plot autocorrelogram, average waveform,

wave_fig=figure('Name', 'Average waveform for each cluster', 'color' ,'w', 'units', 'centimeters', 'pos', [1, 2, 25, 18]);
p = uipanel('Parent',wave_fig,'BorderType','none', 'BackgroundColor', 'w'); 
p.Title = 'Average waveforms'; 
p.TitlePosition = 'centertop'; 
p.FontSize = 12;
p.FontWeight = 'bold';

auto_fig=figure('Name', 'Autocorrelogram for each cluster', 'color', 'w', 'units', 'centimeters', 'pos', [3, 2, 25, 18]);
p2 = uipanel('Parent',auto_fig,'BorderType','none','BackgroundColor', 'w'); 
p2.Title = ['Autocorrelograms for each cluster,' num2str(params.binsize) 'ms bins']; 
p2.TitlePosition = 'centertop'; 
p2.FontSize = 12;
p2.FontWeight = 'bold';


for iU=1:length(GDcells)
    iUnit=GDcells(iU);
    ts_= spikeStruct.timesSorted{iUnit}; 
    n_spk=length(ts_); %how many spikes in the cluster.
    
    %estimate the baseline FR by convolving with a 50ms gaussian (may
    %remove this)
%   unit_sdf = unitSDF(spikeStruct, iUnit, [bl_start, bl_end], 50);
%   avFR=mean(unit_sdf);
   
    % alternatively, estimate firing rate by binning and
    % averaging, using the binsize set in params, at the same time as doing
    % the acor.
    binsize = params.binsize/1000;    
    tbin_edges = min_t:binsize:max_t;  %look at FR properties in baseline only
    tbin_centers = tbin_edges(1:end-1)+binsize/2;
    [bincount, ~] = histcounts(ts_, tbin_edges);  

     %calculate the autocorrelation:
    maxlag=ceil(xcor_window/binsize);   %set the lags we care about (in number of bins)
    [r, lags]=xcorr(bincount, maxlag);   
    
    % calculate mean firing rate
    avFR=mean(bincount);
    
    %plot
    figure(auto_fig)
    sp=subplot(sq_, sq_, iU, 'Parent', p2);

    r(find(~lags))=0; %don't plot the zero autocor as will be a very high spike, distorting the plot
    plot(lags*binsize, r)
    c=find(lags*binsize>-xcor_window & lags*binsize<xcor_window);
    plotmax=1.01 * max(r(c));
    xlim([-xcor_window, xcor_window]);
    ylim([0, plotmax]);
    xlabel('Time (s)')
    ylabel('autocorr')
    title(['Clu #', num2str(iUnit)])
    
    %pull out the average waveform for each cluster, with its standard
    %error.
    wave=spikeStruct.av_waveform{iUnit};
    
    if isfield(spikeStruct, 'std_waveform')
       wfstd=spikeStruct.std_waveform{iUnit};
       if isfield(spikeStruct, 'nWFs_extracted')
          n_spk_wf=spikeStruct.nWFs_extracted(iU); %if we saved the # spikes used to calculate wfs
          wf_sem=wfstd/(n_spk_wf^0.5);
       else
           if n_spk <2000  %tif we don't actually know, assume it used everything it could up to a limit of 2000.
               wf_sem=wfstd/(n_spk^0.5);
           else
               wf_sem=wfstd/(2000^0.5);
           end
       end    
    else
        wf_sem=zeros(1,length(wave)); %don't plot it if it doesnt exist.
    end
    
    %plot the average waveform:
    wave_time=1000/fs *(1:length(wave)); %convert to milliseconds
    figure(wave_fig)  
    thissp=subplot(sq_, sq_, iU, 'Parent', p);
    shadedErrorBar(wave_time, wave, wf_sem, 'b', .5)
    xlabel('ms')
    ylabel('\mu V')
    title({['Clu #', num2str(iUnit)],...
        ['av. FR=', num2str(avFR,2), 'Hz']})
    ax = gca;
    ax.FontSize = 8;
    hold on;       
  
end

end

