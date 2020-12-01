% A plotting script that starts with a completed spikestruct, plots mean of
% multiple laser and footshock events for all clusters identified in order to identify LC cells.
% Anna Sales 2017 Last update June 2018


%% Load in the spikeStruct.
 datapath=pwd; %dont forget backslash
% load([datapath 'spikeStruct.mat']);

% Path to ADC file with pedal on
ADC_fn=[datapath '\100_ADC3.continuous'];
% Path to continuous file with ECG recording.
ECG_fn=[datapath '\100_CH34.continuous'];

%% specify TTLs to plot

pedal_ts=pedalOnOffs(ADC_fn);
TTLchan1=7;  %this is the channel to take TTLs from
fsts=spikeStruct.TTLs.digital{TTLchan1};
% fsts(2:2:end)=[];

inds2=find( diff(fsts)>1.9 & diff(fsts) <2.1)+1;
inds2=[inds2(1)-1; inds2];
fsts_2Hz=fsts(inds2);

inds025=find( diff(fsts)>3.9 & diff(fsts) <4.1)+1;
inds025=[inds025(1)-1; inds025];
fsts_025Hz=fsts(inds025);

nlabels=length(spikeStruct.TTLs.manual.TTL_labels);
labels_to_keep=[2,7,8,10,11,12,14,23,24,25,26,27,30];
labels=cell(1, length(labels_to_keep));
labels(1:length(labels_to_keep))=spikeStruct.TTLs.manual.TTL_labels(labels_to_keep)

labeltimes=spikeStruct.TTLs.manual.TTL_times(labels_to_keep)

%% Pull out useful things for plotting from the spikeStruct:

nclusts=spikeStruct.nclusts;  %number of clusters
fs=spikeStruct.sample_rate;   % sampling rate
newWFs=spikeStruct.allchanWFs;  %waveforms across channels, for each cluster
c_channel=spikeStruct.c_channel;  %centre channel for the cluster
av_waveform=spikeStruct.av_waveform;  %av waveform on the centre channel
plot_pos=spikeStruct.plot_pos;     % depth of each cluster
bl_start=spikeStruct.baseline_st;  %baseline info
bl_end=spikeStruct.baseline_end;
min_t=spikeStruct.timeRange(1);    %time range of recording
max_t=spikeStruct.timeRange(2);
sq_=ceil(nclusts^0.5) ; %for calculating the number of subplots required

%set up some labels for plots below.
for pos=1:1:length(plot_pos)
    unit_test(pos)=find(plot_pos==pos); %The unit that is in the pos-th position on the plot
    chan_=spikeStruct.c_channel(unit_test(pos));   
    tt=['Clu ', int2str(unit_test(pos)), '  c chan= ', int2str(spikeStruct.c_channel(unit_test(pos)))];
    ticklabs{pos}=tt;
end

ticker=11:10:(10*(nclusts+1));

%% Split spikes around specified events into two groups, splitting by latency threshold
%splitting threshold in seconds:
% (Spikes occuring within this time of the event will be in group 1; 
% those occuring after (but within win(2)) will be in group 2)

lat_threshold=0.074; %splitting
min_lat=0.06;  %minimum time for spike
max_lat=0.079;   % max time for spike

clust_to_check=5; %the unit number we are trying to split.
win=[0.25, 0.1] ; %specify a window, time before and after event to consider for pinches
binwin=0.002;
foot_sp_fig      = figure('color','w','NumberTitle','off', 'name','Spiking around footshock TTLs', 'units', 'centimeters', 'pos',[5 2 24 17]);
% foot_spTabGroup = uitabgroup(foot_sp_fig,'TabLocation','Left');
fsts_keep=fsts([290:480]);
TTL_to_plot=fsts_keep;  %update this as needed. These are the TTLs to plot.
nTTL=length(TTL_to_plot);
short_lats=[];
long_lats=[];

for iUnit=clust_to_check; %the unit number we are trying to split.

    
  all_short_lat=[];
  all_long_lat=[];
  
  figure(foot_sp_fig)
%   footunit_tab_sp = uitab(foot_spTabGroup, 'Title', ['Cluster ' num2str(iUnit)],'BackgroundColor','w');
%   axes('Parent',footunit_tab_sp);

  %pull out the relevant spike times
  ts_= spikeStruct.timesSorted{iUnit};

  event_ts=[];
  spk_count_all=[];
  
  for iTTL=1:nTTL %footshocks - clean this up later, this is a long winded way of doing this

      event_ts=TTL_to_plot(iTTL);

      win_st=event_ts-win(1);
      win_end=event_ts+win(2);

      t_ind1=find(ts_>=win_st & ts_<=win_end);
      ts_window=ts_(t_ind1);  %store spike times, relative to event
      
      short_lat_=ts_window(ts_window-event_ts<lat_threshold & ts_window-event_ts>min_lat & ts_window-event_ts<max_lat);
      long_lat_=ts_window(ts_window-event_ts>lat_threshold & ts_window-event_ts>min_lat & ts_window-event_ts<max_lat);
     
      all_short_lat=[all_short_lat; short_lat_];
      all_long_lat=[all_long_lat; long_lat_];
   
      d=subplot(2,1,1);
  
      reps=5;
      if reps==length(ts_window);
          reps=6;  %had to put this in because the plot will mess up if ts_plot is a square!
      end

      ts_plot_short=repmat(short_lat_-event_ts, 1, reps);  %NB if this is a square matrix the plot will mess up as it'll go along wrong dim
      ts_plot_long=repmat(long_lat_-event_ts, 1, reps);  
      y_marks_=linspace(-0.3,0.3, reps) + (iTTL) ; %centres at  1,2,3 etc for each trial  #

      if ts_plot_short
          plot(ts_plot_short, y_marks_', 'Color', rgb('SeaGreen'));         
          hold on
      end
      
      if ts_plot_long
          plot(ts_plot_long, y_marks_', 'b');
      end
      hold on
      xlim([-0.05 0.1]);

  end

   xlabel('Time (s)')
  ylabel('Trial #')
  set(gca, 'FontSize', 11);
  title(['Spikes around footshock, unit #' num2str(clust_to_check)], 'FontWeight', 'normal')    
  yticks(0:10:nTTL);
  xlim([-0.05, 0.1]);
  ylim([0, nTTL+1]);
  plot( zeros(1,2), [0, nTTL+0.5], 'r', 'LineWidth', 1)  %use an odd number in the plotting as laser stamps come in pairs - that way it'll never be square

 
  d.Position=[0.12, 0.1, 0.7, 0.84];
  long_lats=all_long_lat;
  short_lats=all_short_lat;
end

%% Plot waveforms for selected short/long group
raw_data_bin='dataEditedALL.bin';

spike_times_short=short_lats;
spike_times_long=long_lats;

all_spk_ts=spikeStruct.st;
all_spk_IDs=spikeStruct.clu;

%now pull out the cluster ID for each of these spikes (not the same as the
%ID of the clustered unit! This is a template ID generated by kilosort.
short_inds=find(ismember(all_spk_ts, spike_times_short));
long_inds=find(ismember(all_spk_ts, spike_times_long));

% now make sure these times are all associated with the right cluster (as
% rarely we get spikes on multiple clusters at exactly the same time)
short_inds(spikeStruct.clu(short_inds)~=spikeStruct.cids(clust_to_check))=[];
long_inds(spikeStruct.clu(long_inds)~=spikeStruct.cids(clust_to_check))=[];

spike_IDs_short=all_spk_IDs(short_inds);
spike_IDs_long=all_spk_IDs(long_inds);

spike_inds_short=round(short_lats*fs);  %relative to raw continuous data
spike_inds_long=round(long_lats*fs);

wfparams.dataDir=[pwd '\'];; %where is raw data
wfparams.fileName=raw_data_bin;
wfparams.nCh = 32;
wfparams.dataType='int16';
wfparams.wfWin = [-40 41];  %number of samples before and after spike peak
wfparams.nWf = 2000;   %number of waveforms to return (if they are there)
wfparams.nBad=0; %nobad channels.
wfparams.cids=spikeStruct.cids;

fprintf('\n Extracting waveforms for short latency spikes...')
wfparams.spikeTimes=spike_inds_short;  %
wfparams.spikeClusters = spike_IDs_short;  %IDs of each spike, when all spikes are listed in one vector
wf_short=getWaveForms_with_std_and_bad(wfparams);   %extract the waveforms

fprintf('\n Extracting waveforms for long latency spikes...')
wfparams.spikeTimes=spike_inds_long;  %
wfparams.spikeClusters = spike_IDs_long;  %IDs of each spike, when all spikes are listed in one vector
wf_long=getWaveForms_with_std_and_bad(wfparams);   %extract the waveforms

wave_time=1000/fs *(1:size(wf_long.waveFormsMean, 3)); %convert to milliseconds
cent_chan=spikeStruct.c_channel(clust_to_check)

num_waves_long=sum(~isnan(wf_long.spikeTimeKeeps));  %how many waveforms to consider.
num_waves_short=sum(~isnan(wf_short.spikeTimeKeeps));

figure('color', 'w')
short_plot=shadedErrorBarLight(wave_time, wf_short.waveFormsMean(1,cent_chan, :),wf_short.waveFormsSTD(1,cent_chan, :),{'g','markerfacecolor',rgb('SeaGreen')},1)
hold on
long_plot=shadedErrorBarLight(wave_time, wf_long.waveFormsMean(1,cent_chan, :),wf_long.waveFormsSTD(1,cent_chan, :), 'b',1)
xlabel('Time (s)')
ylabel('Voltage \muV')
legend([long_plot.mainLine, short_plot.mainLine], {'Long latency', 'Short latency'})
%% Plot across chans


dist=100; %channel spread aorund centre channel to plot

%SHORTS
%Get waveforms out of wf_short
wfs_short=squeeze(squeeze(wf_short.waveFormsMean));
wfs_short_sem=squeeze(wf_short.waveFormsSTD);
wfs_short_n_extracted=length(find(~isnan(wf_short.spikeTimeKeeps)));
wave_time=1000/spikeStruct.sample_rate *(1:size(wfs_short,2)); %convert to milliseconds

%LONGS
%Get waveforms out of wf_longs
wfs_long=squeeze(squeeze(wf_long.waveFormsMean));
wfs_long_sem=squeeze(wf_long.waveFormsSTD);
wfs_long_n_extracted=length(find(~isnan(wf_long.spikeTimeKeeps)));

cent_chan=spikeStruct.c_channel(clust_to_check);  %the centre channel for this unit

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
pltcol_short=rgb('Green');
pltcol_long=rgb('Blue')
for sb=1:length(chans_within_range)
    chan_=chans_within_range(sb);
       
    %put the plot in the correct place.
    chanPlots(sb)=subplot('Position', [plotPos(sb,1), plotPos(sb,2), 0.18, 0.1]);
    long_plot=shadedErrorBarLight(wave_time, wfs_long(chan_, :), wfs_long_sem(chan_,:)./sqrt(wfs_long_n_extracted), 'b', 1);
    hold on
    short_plot=shadedErrorBarLight(wave_time, wfs_short(chan_, :), wfs_short_sem(chan_,:)./sqrt(wfs_short_n_extracted), 'g', 1);
    short_plot.mainLine.Color=rgb('LimeGreen');
    short_plot.patch.FaceColor=rgb('LimeGreen');
    hold on
    if chan_==cent_chan  
        text(2,100, 'Centre channel', 'FontSize', 7)
    end
    
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
text(3.3, 150, ['Unit: ' num2str(clust_to_check)]);
text(3.3, 0, ['Green=short lat., blue=long lat.']);
text(3.3, -80, ['Showing channels within ' num2str(dist) '\mum']);