% A plotting script that starts with a completed spikestruct, plots mean of
% multiple laser and footshock events for all clusters identified in order to identify LC cells.
% Anna Sales 2017 Last update June 2018


%% Load in the spikeStruct.
 datapath='D:\NTS_Pabitra\120820\Rec2\'; %dont forget backslash
 load([datapath 'spikeStruct.mat']);
 
 ADC_fn='D:\NTS_Pabitra\120820\Rec2\ADC\100_ADC2.continuous';
%% Plot out where all the clusters are in space on the probe
probe_fig = plotProbe(spikeStruct)

%% Plot waveforms and autocorrelograms

%specify parameters
pms.cell_list=1:spikeStruct.nclusts;
pms.binsize=2;  %binsize, in ms
pms.window=250;  %window size, in ms.
pms.baseline = 1 ;

[wavefig, autofig] = plot_waveforms_acors(spikeStruct, pms)

%% Plot xcors, to check isolation - NB don't plot too many as it's slow
pms2.cell_list=[1,2,3]%spikeStruct.nclusts];
pms2.binsize=5;  %binsize, in ms
pms2.window=500;  %window size, in ms.
pms2.baseline = 1 ;

[corfig] = cell_xcors(spikeStruct, pms2)

%% specify TTLs to plot

TTLchan=6;  %this is the channel to take TTLs from - 6 for footshock
TTLts=spikeStruct.TTLs.digital{TTLchan};
TTLts(2:2:end)=[]; %Because we have an 'on' and an 'off' and usually only want 'on'


%% Test out the laser widget
laserchan=5;  %this is the channel to take TTLs from
laserts=spikeStruct.TTLs.digital{laserchan};
if laserts
    [stims, stims_by_type]=laserTTLwidget(laserts);
end
%% plot waveforms during laser stim and compare to those from other times

short_pulses=stims_by_type{1};
pulse_times=short_pulses(:,1);  %times of pulses
duration=5; %pulse duration in ms
wf_fig=wf_during_laser(spikeStruct,pulse_times, duration, 'dataALL.bin' );

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
    ticklabs_clu{pos}=tt;
end

ticker=11:10:(10*(nclusts+1));

%% Plotting activity around a TTL - change fields as required
  
win=[0.25, 0.5] ; %specify a window, time before and after event to consider for shocks
binwin=0.01;
foot_sp_fig  = figure('color','w','NumberTitle','off', 'name','Spiking around footshock TTLs', 'units', 'centimeters', 'pos',[5 2 24 17]);
foot_spTabGroup = uitabgroup(foot_sp_fig,'TabLocation','Left');

TTL_to_plot=TTLts;  %update this as needed. These are the TTLs to plot.
%will need amending if we want to plot laser stuff.
nfs=length(TTL_to_plot);

for iUnit=1:nclusts

  figure(foot_sp_fig)
  footunit_tab_sp = uitab(foot_spTabGroup, 'Title', ['Cluster ' num2str(iUnit)],'BackgroundColor','w');
  axes('Parent',footunit_tab_sp);

  %pull out the relevant spike times
  ts_= spikeStruct.timesSorted{iUnit};

  event_ts=[];
  spk_count_all=[];
  for iTTL=1:nfs %footshocks - clean this up later, this is a long winded way of doing this

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
      xlim([-0.5 0.5]);
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
      xlim([-0.5 0.5]);

  end

  xlim([-0.5 0.5])
  xlabel('')
  ylabel('Trial #')
  set(gca, 'FontSize', 11);
  title('Spikes around footshock', 'FontWeight', 'normal')    
  yticks(1:nfs);
  xlim([-win(1), win(2)]);
  ylim([0, nfs+1]);
  plot( zeros(1,2), [0, nfs+0.5], 'r', 'LineWidth', 1)  %use an odd number in the plotting as laser stamps come in pairs - that way it'll never be square

  %get the means / sums

  %plot mean firing rate
  binned_mean_fr=mean(spk_count_all,1);
  f=subplot(2,1,2);
  bar(t_plot,binned_mean_fr);
  hold on
  xlabel('Time (s)')
  ylabel('Mean spike count')
  text(0.3, 2,[num2str(1000*binwin) 'ms bins'])
  xlim([-win(1), win(2)]);
  maxploty=max(binned_mean_fr)+1;
  minploty=0;
  plot( zeros(1,7), linspace(0, maxploty, 7), 'r', 'LineWidth', 1)  %use an odd number in the plotting as laser stamps come in pairs - that way it'll never be square
  ylim([0,2.5])
  set(gca, 'FontSize', 11);


  d.Position=[0.12, 0.35, 0.7, 0.6];
  f.Position=[0.12, 0.075, 0.7, 0.2];
end

%% Smooth firing rate repsonses on the probe by depth, with TTLs. 
cellList=1:nclusts;
fs=spikeStruct.sample_rate;  %sampling rate
c_channel=spikeStruct.c_channel;  %centre channel for the cluster
min_t=spikeStruct.timeRange(1);  %min, max spikes in the recording
max_t=spikeStruct.timeRange(2);
bl_start=spikeStruct.baseline_st;  %baseline times, if needed.
bl_end=spikeStruct.baseline_end;

c_channel2=c_channel(cellList);
[vals, indys]=sort(c_channel2);

%% FR plot - one for each cell, plotted in correct depth order.

timePlot=[0 30]; % time to plot, in seconds.
binsize=0.2; %in seconds
ylimits=[0, 33];
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
ticklabs_clu=cell(1,nclusts2); %plot labels for clusters
ticklabs_ch=cell(1,nclusts2);   %plot labels for channels
%set up some labels for plots below.
for pos=1:1:length(plot_pos2)
    unit_test(pos)=find(plot_pos2==pos); %The unit that is in the pos-th position on the plot
    chan_=c_channel2(unit_test(pos));   
    tt=['Clu. ' int2str(cellList(unit_test(pos))) ' Ch. ' int2str(c_channel2(unit_test(pos)))];
    ticklabs_clu{pos}=[tt];
    ticklabs_ch{pos}=int2str(c_channel2(unit_test(pos)));
end 


%check if there are laser TTLs included in this period - and mark on the
%plot if there are.

TTL_inc=TTLts(TTLts>timePlot(1)& TTLts<timePlot(2));

fr_fig= figure('color','w','NumberTitle','off', 'name',' Unit firing', 'units', 'centimeters', 'pos',[5 2 24 18.5], 'Color', 'white');
p6 = uipanel('Parent',fr_fig,'BorderType','none'); 

figure(fr_fig)
ADC_plot=subplot('Position',[0.064 0.075 0.85 0.15],  'Parent', p6);

FRsubplots={};
sp_height=0.61/length(cellList);
if sp_height > 0.3
    sp_height=0.24
end


for u=1:length(cellList)
    
    FRsubplots{u}=subplot('Position',[0.064 0.26+((plot_pos2(u)-1)*1.2*sp_height) 0.85 sp_height],  'Parent', p6);
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
    spk_count = spk_count(1:end-1)/binsize;
    %plot
    mybar=bar(t_plot, spk_count);
    mybar.FaceColor=[0 0 1];
    mybar.FaceAlpha=0.7;
    mybar.LineWidth=1;
    mybar.EdgeColor='none';
    
    xlim([timePlot(1), timePlot(2)]);
    ylim(ylimits);
    hold on
 
    left_label=ylabel('' );
    hold on
    
    if plot_pos2(u) >1
        set(gca,'Xticklabel',[]) ;
        xlabel('');
    end
    
    if plot_pos2(u) == floor(nclusts/2)
        text_lab=text(-0.05,1,['Spike count / binsize (' num2str(binsize) 'ms)'],'Units', 'Normalized', 'FontSize', 10);
        text_lab.Rotation=90;
    end
       
    if mod(plot_pos2(u),2)
         set(gca,'Yticklabel',[]) ;
    end
    aa=gca;
    aa.FontSize=8;
    ymax=ylimits(2);
    ymin=0;
    left_label.Color='k';
    left_label.Position=[198.5, 3, -1];
    left_label.FontSize=8;

    if ~isempty(TTL_inc)     
        plot(repmat(TTL_inc,1,3), [ymin, ymax, ymax], '-r', 'Color', [1 0 0 0.6], 'LineWidth', 2);
    end
    
    yyaxis right
    label_right=ylabel(ticklabs_clu{plot_pos2(u)}, 'FontSize', 7);
    set(gca,'XTick',[], 'YTick', []);
    set(gca,'YColor','k');
    label_right.Rotation=0;
    label_right.Position=[timePlot(2)+(0.01*range(timePlot)), 1, -1];
    label_right.HorizontalAlignment='left';
end

%
%now plot ECG
[raw_sig, ADC_ts_raw,  ADC_info] = load_open_ephys_data(ADC_fn);  
ADC_tbase=ADC_ts_raw-ADC_ts_raw(1); % timebase should start at zero 

keep_=find(ADC_tbase>timePlot(1)-10 & ADC_tbase<timePlot(2)+10);

ADC_tbase_win=ADC_tbase(keep_);
raw_sig_win=raw_sig(keep_);

% Now do some filtering / editing
zsig=zscore(raw_sig_win);
mn_sg=mean(zsig);
std_sig=std(zsig);
nstd=2;
high_ind=find(abs(zsig)>nstd*std_sig); %set anything above 2stds back to the mean (remove high ECG readings)
zsig(high_ind)=0;
%take the moving mean, to loose higher frequency ECG
resps=movmean(abs(zsig),8000);

%now high pass filter to move slow drifts in the baseline
HPlim=0.2;   %filter to take out slow drifts
[Db,Da]=butter(2,HPlim/(0.5*fs), 'high' );
HP_resps=filtfilt(Db,Da,resps);  %filtered delta

subplot(ADC_plot)
 hold on
 plot(ADC_tbase_win, (raw_sig_win));
 ylim([-4,8]);
 ylabel('ECG');
 yyaxis right;
plot(ADC_tbase_win,HP_resps+0.2, 'r');
 ylim([-0.2,0.5]);
ylabel('Resp');
 xlim([timePlot(1), timePlot(2)]);
 xlabel('Time (s)');
 set(gca,'YColor','r');
%now overlay laser TTL markers, if included in the period specified.
plt=[];
ymin=-5;

if ~isempty(TTL_inc)
    plot(repmat(TTL_inc,1,3), [ymin, ymax, ymax], '-b', 'Color', [1 0 0 0.6], 'LineWidth', 2);
end


xlim([timePlot(1), timePlot(2)])
% title('ECG/resps', 'FontWeight', 'normal', 'FontSize', 10)

%% Now calculate smooth spike density function for individual/all units, and look at 
%  xcor with respiratory cycle.
timeWin=[0, 300];

ADCfs=ADC_info.header.sampleRate;
%downsample to 1000, to match up with the binned spiking (1ms bins)

ADCds=ADCfs/1000; %factor to downsample by
ADC_tbase_ds=downsample(ADC_tbase, ADCds);
ADC_sig_ds=downsample(raw_sig,ADCds);
ADCfs=ADCfs/ADCds;

keep_=find(ADC_tbase_ds>timeWin(1)& ADC_tbase_ds<timeWin(2));
ADC_tbase_win=ADC_tbase_ds(keep_);
raw_sig_win=ADC_sig_ds(keep_);
% Now do some filtering
zsig=zscore(raw_sig_win);
mn_sg=mean(zsig);
std_sig=std(zsig);
nstd=2;
high_ind=find(abs(zsig)>nstd*std_sig);
zsig(high_ind)=0;
resps_win=movmean(abs(zsig),8000);

HPlim=0.2;   %filter to take out slow drifts
[Db,Da]=butter(2,HPlim/(0.5*ADCfs), 'high' );
HP_resps_win=filtfilt(Db,Da,resps_win);  %filtered delta

% quick sanity check
% figure
% plot(ADC_tbase_win, resps_win-mean(resps_win))
% hold on
% plot(ADC_tbase_win, HP_resps_win)
% title('HP filtered resp signal')

% 
figure('Color', 'w', 'Units', 'centimeters', 'Position', [2 2 22 19])
nplots=ceil(sqrt(nclusts));
sdf_store=[];
pl=1;
max_lag=550; %max lag to consider, in ms
ylimxcor=[-0.1 0.1];
for h=1:nclusts

     [sdf, sdf_tbase] = unitSDF(spikeStruct, h, timeWin ,100);
     sdf_store(h,:)=sdf;
    %NB xcor is between fixed first vector, and lagged version of second.
    % So lag is relative to SECOND signal (here, sdf)

    if sum(sdf)
        [xc, lags]=xcorr(HP_resps_win,sdf, 'coeff', max_lag); %'lags' will run to/from +/- (centres(end)-1)
        subplot(nplots, nplots, pl)
        plot(lags, xc)
        hold on
        plot([0 0], [-1 1], ':k');
        plot([-max_lag, max_lag], [0 0], ':k');
        aa=gca;
        aa.FontSize=7;
        if h==nclusts
           xl=xlabel('Time (ms)', 'FontSize', 7);
        end
        yl=ylabel('XC', 'FontSize', 8);
        yl.Position(1)=-1.3*max_lag;
        title(['Clust #' num2str(h)], 'FontWeight', 'normal', 'FontSize', 9)
        ylim(ylimxcor)
        xlim([-max_lag, max_lag])
        pl=pl+1;
    end
    
    
end
text(2*max_lag, 0, {'X-cor between respiratory signal' 'and spike density over 5m baseline'})

%% Sanity check

figure
subplot(2,1,1)
plot(ADC_tbase_win, resps-mean(resps))
hold on
plot(ADC_tbase_win, HP_resps)
title('HP filtered resp signal')
subplot(2,1,2)
plot(sdf_tbase, sdf_store(1,:))
title('Example spike SDF')
%% Now merge all the spike trains, see what happens.

[sdf_all, sdf_tbase_all] = unitSDF_all(spikeStruct, 1:nclusts, timeWin ,100);

 max_lag=200; %max lag to consider, in ms
[xc, lags]=xcorr(HP_resps,sdf_all,  'coeff', max_lag); %'lags' will run to/from +/- (centres(end)-1)
figure
plot(lags, xc)
hold on
plot([0 0], [0 1], ':k')
plot([-max_lag, max_lag], [0 0], ':k')
xlabel('Time (ms)', 'FontSize', 8)
ylabel('XC', 'FontSize', 8)
title('Merged spike train', 'FontWeight', 'normal', 'FontSize', 8)
ylim([0, 0.3])


    