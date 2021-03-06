% A plotting script that starts with a completed spikestruct, plots mean of
% multiple laser and footshock events for all clusters identified in order to identify LC cells.
% Anna Sales 2017 Last update June 2018

addpath(genpath('D:\Code\MATLAB\AnalysingPHYoutput\GENERIC'));
%% Load in the spikeStruct.
%  datapath='D:\Versus\2807\2807oepformat\20200728\Rec1\CommonAvRef\'; %dont forget backslash
%  load([datapath 'spikeStruct.mat']);
%% Plot out where all the clusters are in space on the probe
probe_fig = plotProbe(spikeStruct)

%% Plot waveforms and autocorrelograms

%specify parameters
pms.cell_list=[1:spikeStruct.nclusts];
pms.binsize=1;  %binsize, in ms
pms.window=150;  %window size, in ms.
pms.baseline = 1 ;

[wavefig, autofig] = plot_waveforms_acors(spikeStruct, pms)

%% Plot across multiple chans
for u=1:spikeStruct.nclusts
    multichan_wfs=plot_across_chans(spikeStruct, u) 
    savefig( multichan_wfs, [pwd '\unit' num2str(u) 'wf_across_chans.fig'] );  
end

 
%% Plot xcors, to check isolation - NB don't plot too many as it's slow
pms2.cell_list=1:spikeStruct.nclusts;
pms2.binsize=1;  %binsize, in ms
pms2.window=150;  %window size, in ms.
pms2.baseline = [] ; %leave empty if no baseline!

[corfig] = cell_xcors(spikeStruct, pms2)

%% specify TTLs to plot

TTLchan1=3;  %this is the channel to take TTLs from
fsts=spikeStruct.TTLs.digital{TTLchan1};
fsts(2:2:end)=[];
% fsts(fsts<270)=[];
% 
% fsts_2Hz=fsts
% fsts_2Hz(diff(fsts_2Hz)>0.55)=[];
% fsts_2Hz(end)=[];

% fsts_025Hz=fsts
% fsts_025Hz(diff(fsts_025Hz)<4| diff(fsts_025Hz) >4.5 )=[];


TTLchan2=7;
pedal_ts=spikeStruct.TTLs.digital{TTLchan2};
sep_events_inds=find(diff(pedal_ts)>0.2); %separate pedal pushes
ped_ev=[pedal_ts(1);pedal_ts(sep_events_inds+1)]
ped_starts=ped_ev(1:2:end);
ped_ends=ped_ev(2:2:end);
% mants(2:2:end)=[];

% labels_to_keep=[10,15,16,18,19,20,22,24,28];
% labels_to_keep=[3,4,5,6,7,13,14,15,16,17,18,20,22]
labels_to_keep=[1,2,3,4,5]
labels=cell(1, length(labels_to_keep));
labels(1:length(labels_to_keep))=spikeStruct.TTLs.manual.TTL_labels(labels_to_keep)

labeltimes=spikeStruct.TTLs.manual.TTL_times(labels_to_keep)


% Need to edit so we only consider the first one if it's not a laser
% % event!!
% TTLset1=TTLts(find(TTLts<labeltimes(3)));
% TTLset2=TTLts(find(TTLts>labeltimes(3) & TTLts < labeltimes(4)));
% TTLset3=TTLts(find(TTLts > labeltimes(4)));

%% Test out the laser widget
% laserchan=3;  %this is the channel to take TTLs from
% laserts=spikeStruct.TTLs.digital{laserchan};
% [stims, ~]=laserTTLwidget(laserts);

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

%% Plotting activity around a TTL - change fields as required
  
win=[0.25, 0.5] ; %specify a window, time before and after event to consider for pinches
binwin=0.01;
foot_sp_fig      = figure('color','w','NumberTitle','off', 'name','Spiking around footshock / pinch TTLs', 'units', 'centimeters', 'pos',[5 2 24 17]);
foot_spTabGroup = uitabgroup(foot_sp_fig,'TabLocation','Left');

TTL_to_plot=fsts%fsts_2Hz;  %update this as needed. These are the TTLs to plot.
nTTL=length(TTL_to_plot);

for iUnit=1:nclusts

  figure(foot_sp_fig)
  footunit_tab_sp = uitab(foot_spTabGroup, 'Title', ['Cluster ' num2str(iUnit)],'BackgroundColor','w');
  axes('Parent',footunit_tab_sp);

  %pull out the relevant spike times
  ts_= spikeStruct.timesSorted{iUnit};

  event_ts=[];
  spk_count_all=[];
  for iTTL=1:nTTL %footshocks - clean this up later, this is a long winded way of doing this

      event_ts=TTL_to_plot(iTTL);

      win_st=event_ts-win(1);
      win_end=event_ts+win(2);

      t_ind1=find(ts_>=win_st & ts_<=win_end);
      ts_window=ts_(t_ind1);  %store all the data that's been cut.

      % convolve the spike train in the window with Gaussian kernel to
      % estimate FR
      tbin_edges = win_st:binwin:win_end;

      if iTTL==1 %store a time vector for plotting
          tbin_centers = tbin_edges(1:end-1)+binwin/2;
          t_plot=tbin_centers-event_ts;
      end

      %bin and histogram
      spk_count = histc(ts_,tbin_edges);
      spk_count = spk_count(1:end-1);
      spk_count_all(iTTL,:)=spk_count;

      d=subplot(2,1,1);
      xlim([-0.05 0.1])
      reps=5;
      if reps==length(ts_window);
          reps=6;  %had to put this in because the plot will mess up if ts_plot is a square!
      end

      ts_plot=repmat(ts_window-event_ts, 1, reps);  %NB if this is a square matrix the plot will mess up as it'll go along wrong dim
      y_marks_=linspace(-0.3,0.3, reps) + (iTTL) ; %centres at  1,2,3 etc for each trial  #

      if ts_window
          plot(ts_plot, y_marks_', 'k');
      end
      hold on
      xlim([-0.05 0.1]);

  end

 
  xlabel('')
  ylabel('Trial #')
  set(gca, 'FontSize', 11);
  title('Spikes around footshock', 'FontWeight', 'normal')    
  yticks(0:10:nTTL);
  xlim([-0.05, 0.1]);
  ylim([0, nTTL+1]);
  plot( zeros(1,2), [0, nTTL+0.5], 'r', 'LineWidth', 1)  %use an odd number in the plotting as laser stamps come in pairs - that way it'll never be square

  %get the means / sums

  %plot mean firing rate
  binned_mean_fr=mean(spk_count_all,1);
  f=subplot(2,1,2); 
  xlabel('Time (s)')
  ylabel('Mean spike count')
  if binned_mean_fr
      bar(t_plot,binned_mean_fr);
      hold on
    
      text(0, 0.8,[num2str(1000*binwin) 'ms bins'])
      xlim([-win(1), win(2)]);
      maxploty=max(binned_mean_fr)+1;
      minploty=0;
      plot( zeros(1,7), linspace(0, maxploty, 7), 'r', 'LineWidth', 1)  %use an odd number in the plotting as laser stamps come in pairs - that way it'll never be square
      ylim([0,1])
      set(gca, 'FontSize', 11);
      xlim([-0.05 0.1])
  end
  d.Position=[0.12, 0.35, 0.7, 0.6];
  f.Position=[0.12, 0.075, 0.7, 0.2];
end

%% Smooth firing rate repsonses on the probe by depth, with TTLs. 

cellList=[1,2,3]
timePlot=[0,440]; % time to plot, in seconds.
binsize=1; %in seconds

fs=spikeStruct.sample_rate;  %sampling rate
c_channel=spikeStruct.c_channel;  %centre channel for the cluster


c_channel2=c_channel(cellList);
[vals, indys]=sort(c_channel2);

%work out plotting positions so that the cells are plotted as they appear
%on the probe, i.e. cells detected ventrally on the probe are at the bottom
%of the window.

plot_pos2=[];
for iUnit = 1:1:length(cellList)
    
    [~,plot_pos2(iUnit)]=find(indys==iUnit) ; % a plot position for each cluster, with lower (deeper) channels getting lower numbers
    
    %plot_pos is a vector n_clusts long with a position for each unit,
    %provided in the same order as usual.
end

nclusts=length(cellList);
ticker=11:10:(10*(nclusts+1));
ticklabs={nclusts};

%set up some labels for plots below.
for pos=1:1:length(plot_pos2)
    unit_test(pos)=find(plot_pos2==pos); %The unit that is in the pos-th position on the plot
    chan_=c_channel2(unit_test(pos));   
    tt=['Clu ', int2str(cellList(unit_test(pos))), ' ch ', int2str(c_channel2(unit_test(pos)))];
    ticklabs{pos}=tt;
end 


% FR plot - one for each cell, plotted in correct depth order.

c_channel2=c_channel(cellList);
[vals, indys]=sort(c_channel2);

%work out plotting positions so that the cells are plotted as they appear
%on the probe, i.e. cells detected ventrally on the probe are at the bottom
%of the window.

plot_pos2=[];
for iUnit = 1:length(cellList)    
    [~,plot_pos2(iUnit)]=find(indys==iUnit) ; % a plot position for each cluster, with lower (deeper) channels getting lower numbers    
    %plot_pos is a vector n_clusts long with a position for each unit,
    %provided in the same order as usual.
end

nclusts2=length(cellList);
ticker=11:10:(10*(nclusts2+1));
ticklabs=cell(1,nclusts2);

%set up some labels for plots below.
for pos=1:1:length(plot_pos2)
    unit_test(pos)=find(plot_pos2==pos); %The unit that is in the pos-th position on the plot
    chan_=c_channel2(unit_test(pos));   
    tt=[{['Clu ', int2str(cellList(unit_test(pos)))]; [' ch ', int2str(c_channel2(unit_test(pos)))]}];
    ticklabs{pos}=tt;
end 


%check if there are laser TTLs included in this period - and mark on the
%plot if there are.

fsts_inc=fsts(fsts>timePlot(1)& fsts<timePlot(2));
fs2Hz_inc=fsts_2Hz(fsts_2Hz>timePlot(1)& fsts_2Hz<timePlot(2));
fs025Hz_inc=fsts_025Hz(fsts_025Hz>timePlot(1)& fsts_025Hz<timePlot(2));

ped_table=[ped_starts, ped_ends]
ped_inc=ped_ev(ped_ev>timePlot(1)& ped_ev<timePlot(2))

fr_fig= figure('color','w','NumberTitle','off', 'name',' Unit firing', 'units', 'centimeters', 'pos',[5 3 24 15], 'Color', 'white');
p6 = uipanel('Parent',fr_fig,'BorderType','none'); 

figure(fr_fig)

FRsubplots={};
sp_height=0.6/length(cellList);
if sp_height > 0.3
    sp_height=0.24
end


for u=1:length(cellList)
    
    FRsubplots{u}=subplot('Position',[0.065 0.2+((plot_pos2(u)-1)*1.1*sp_height) 0.88 sp_height],  'Parent', p6);
    
    iUnit=cellList(u);
    ts_= spikeStruct.timesSorted{iUnit}; 
    ts_ind=find(ts_>=timePlot(1)& ts_<=timePlot(2));
    ts_window=ts_(ts_ind);  %store all the data that's been cut.
    tbin_edges = timePlot(1):binsize:timePlot(2);

    if u==1 %store a time vector for plotting
        tbin_centers = tbin_edges(1:end-1)+binsize/2;
        t_plot=tbin_centers;
    end
    
    %bin and histogram
    spk_count = histc(ts_,tbin_edges);
    spk_count = spk_count(1:end-1) ./binsize;
    %plot
    mybar=bar(t_plot, spk_count);
    mybar.FaceColor=[0 0 1];
    mybar.FaceAlpha=0.7;
    mybar.LineWidth=1;
    
    
    xlim([timePlot(1), timePlot(2)])
    hold on
    if plot_pos2(u)==1
        xlabel('Time (s)')
    end
    
    if plot_pos2(u)==length(cellList)
        aa=gca;
        text(-20,aa.YLim(2)+10, 'Spk rate /s')
        text(0.8*timePlot(2),aa.YLim(2), 'Green: pedal, red: footshock.')
    end
    ylabel(ticklabs{plot_pos2(u)}, 'FontSize', 8 )
    hold on
    
    if plot_pos2(u) >1
        set(gca,'Xticklabel',[]) 
        xlabel('')
    end
    
    aa=gca;
    ymax=aa.YLim(2)
    ymin=aa.YLim(1)
    
    if fsts_inc
       plot(repmat(fsts_inc,1,3), [ymin, ymax, ymax], '-r', 'Color', [1 0 0 0.1], 'LineWidth', 1);
    end
%     if fs2Hz_inc
%        plot(repmat(fs2Hz_inc,1,3), [ymin, ymax, ymax], '-r', 'Color', [1 0 0 0.1], 'LineWidth', 1);
%     end
%     
%     if fs025Hz_inc
%        plot(repmat(fs025Hz_inc,1,3), [ymin, ymax, ymax], '-r', 'Color', [1 0.5 0.5 0.1], 'LineWidth', 1);
%     end
   
    if ped_inc
       plot(repmat(ped_inc,1,3), [ymin, ymax, ymax], '-g', 'Color', [0 0.8 0.1 0.1], 'LineWidth', 1);
    end
    
end

  MARKERsubplots=subplot('Position',[0.065 0.04 0.88 0.05],  'Parent', p6);
  for r=1:length(labeltimes)
      text(labeltimes(r), 0.1, labels{r}, 'Rotation', 90 , 'Fontsize', 6)
  end
  set(gca,'xtick',[]) 
  set(gca,'ytick',[]) 
  XAxis. TickLength = [0 0];
  xlim([timePlot(1), timePlot(2)])


