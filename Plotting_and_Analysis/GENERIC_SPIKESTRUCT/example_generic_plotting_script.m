% A plotting script that starts with a completed spikestruct, plots mean of
% multiple laser and footshock events for all clusters identified in order to identify LC cells.
% Anna Sales 2017 Last update June 2018


%% Load in the spikeStruct.
 datapath='D:\Versus\2807\2807oepformat\20200728\Rec1\CommonAvRef\'; %dont forget backslash
 load([datapath 'spikeStruct.mat']);

%% Plot out where all the clusters are in space on the probe
probe_fig = plotProbe(spikeStruct)

%% Plot waveforms and autocorrelograms

%specify parameters
pms.cell_list=1:spikeStruct.nclusts;
pms.binsize=5;  %binsize, in ms
pms.window=400;  %window size, in ms.
pms.baseline = 1 ;

[wavefig, autofig] = plot_waveforms_acors(spikeStruct, pms)

%% Plot xcors, to check isolation - NB don't plot too many as it's slow
pms2.cell_list=1:3;
pms2.binsize=5;  %binsize, in ms
pms2.window=500;  %window size, in ms.
pms2.baseline = 1 ;

[corfig] = cell_xcors(spikeStruct, pms2)

%% specify TTLs to plot

TTLchan=5;  %this is the channel to take TTLs from
TTLts=spikeStruct.TTLs.digital{TTLchan};
TTLts(2:2:end)=[]; %Because we have an 'on' and an 'off' and usually only want 'on'
%% Test out the laser widget
laserchan=5;  %this is the channel to take TTLs from
laserts=spikeStruct.TTLs.digital{laserchan};
[stims, stimtypes, stim_utils]=laserTTLwidget(laserts, 0.5);
%% Example - get waveforms during specified times NB you need the raw binary for this.

wfparams.dataDir=[pwd '\'];; %where is raw data
wfparams.fileName='dataALL.bin';
wfparams.nCh = 32;
wfparams.nBad = nChansBad;
wfparams.dataType='int16';
wfparams.wfWin = [-40 41];  %number of samples before and after spike peak
wfparams.nWf = nWF;   %number of waveforms to return (if they are there)

%make a table of times to consider,e.g during TTLs
spikes_duringTTL=[];
clusterIDs_duringTTL=[];  
all_spk_ts=spikeStruct.st;
all_spk_IDs=spikeStruct.clu;

for p=1:length(TTLts) 
    inds_duringTTL=find(all_spk_ts>TTLts(p) & all_spk_ts<TTLts(p)+0.005); %5ms around TTL
    spikes_duringTTL=[spikes_duringTTL,round(all_spk_ts(inds_duringTTL) * fs)']; %spike times of interest, in samples.
    clusterIDs_duringTTL=[clusterIDs_duringTTL, all_spk_IDs(inds_duringTTL)'];
end

wfparams.spikeTimes=spikes_duringTTL';  %
wfparams.spikeClusters = clusterIDs_duringTTL';  %IDs of each spike, when all spikes are listed in one vector

wf=getWaveForms_with_std_and_bad(wfparams);


 %% Plot during specified events, compared to during entire rec


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

TTL_to_plot=fs_set1;  %update this as needed. These are the TTLs to plot.
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
%% Plotting activity around laser
win=[0.5, 0.5] ; %specify a window, time before and after event to consider
 
laser_fr_fig      = figure('color','w','NumberTitle','off', 'name','FR around laser TTLs (gauss kernel)', 'units', 'centimeters', 'pos',[5 2 33 18]);
laser_frTabGroup = uitabgroup(laser_fr_fig,'TabLocation','Left');

laser_sp_fig     = figure('color','w','NumberTitle','off', 'name','Spiking activity around laser TTLs per trial', 'units', 'centimeters', 'pos',[5 2 33 18]);
laser_spTabGroup = uitabgroup(laser_sp_fig,'TabLocation','Left');

sq_=ceil(nclusts^0.5) ; %for calculating the number of subplots required
sq2_=ceil(size(stimtypes,1)^0.5);
iUnit=1;
stimMarkers=stim_utils.Markers; %markers for laser events on plot
laser_labels=stim_utils.labels;  % labels for each laser eventcclc
for iUnit=1:nclusts
        
         %plot each unit on a new tab
         figure(laser_fr_fig)
         unit_tab = uitab(laser_frTabGroup, 'Title', ['Cluster ' num2str(iUnit)],'BackgroundColor','w');     
          axes('Parent',unit_tab);
         
         figure(laser_sp_fig)
         unit_tab_sp = uitab(laser_spTabGroup, 'Title', ['Cluster ' num2str(iUnit)],'BackgroundColor','w');     
          axes('Parent',unit_tab_sp);
         ts_=[];
         
         %pull out the relevant spike times
         ts_= spikeStruct.timesSorted{iUnit};
 
        for stype=1:size(stimtypes,1)
                       
            figure(laser_sp_fig)
            axes('Parent',unit_tab_sp);
            c=subplot(sq2_, sq2_, stype);
            title(laser_labels{stype})      
            stims_=stimtypes{stype};
            nstims=size(stims_, 1);
            gau_sdf_store=[];
            for iTTL=1:nstims
                
                event_ts=stims_(iTTL, 1);   

                win_st=event_ts-win(1);
                win_end=event_ts+win(2);
                 
                ts_ind=find(ts_>=win_st & ts_<=win_end);
                ts_window=ts_(ts_ind);  %store all the data that's been cut.
                                       
                % convolve the spike train in the window with Gaussian kernel to
                % estimate FR
                binsize = 0.001; % select a small bin size for good time resolution
                tbin_edges = win_st:binsize:win_end;
                
                if iTTL==1 %store a time vector for plotting
                    tbin_centers = tbin_edges(1:end-1)+binsize/2;
                    t_plot=tbin_centers-event_ts;
                end
                
                %bin and histogram
                spk_count = histc(ts_,tbin_edges);
                spk_count = spk_count(1:end-1);
                               
                gauss_window = 1./binsize; % 1 second window
                gauss_SD = 0.05./binsize; % 0.05 seconds (50ms) SD
                gk = gausskernel(gauss_window,gauss_SD); 
                gk = gk./binsize; % normalize by binsize
                gau_sdf = conv2(spk_count,gk,'same'); % do the convolution
                gau_sdf_store(iTTL,:)=gau_sdf;
               
                
                %plot the spiking activity on the spike plot
                subplot(c)
                reps=5;
                if reps==length(ts_window);
                    reps=6;  %had to put this in because the plot will mess up if ts_plot is a square!
                end
                 
                 ts_plot=repmat(ts_window-event_ts, 1, reps);  %NB if this is a square matrix the plot will mess up as it'll go along wrong dim  
                 y_marks_=linspace(-0.3,0.3, reps) + (iTTL) ; %centres at  1,2,3 etc for each trial  #
                 if ts_window %only try to plot if there are actually spikes in the window
                    plot(ts_plot, y_marks_', 'k');
                 end
                 hold on            
                         
            end
            
            %add labels / laser markers to the subplot etc.
            subplot(c)
            title(laser_labels{stype})      
            xlabel('Time (s)')
            ylabel('Trial #')
            set(gca,'FontSize', 14);
            yticks(1:nstims);          
            markers=stimMarkers{stype};
            %use an odd number for 'reps' here as should always have an
            %even number of laser markers (on/off TTLs)
            maxploty=nstims+1;
            plot( repmat(markers, 1,7), linspace(0, maxploty, 7), 'Color', [1 0 0 0.25], 'LineWidth', 1)  %use an odd number in the plotting as laser stamps come in pairs - that way it'll never be square
            
            for hh=1:length(markers)/2 %give the laser regions a little shading
                x_patch= [markers((2*hh)-1), markers(2*hh), markers(2*hh), markers(2*hh-1)];
                y_patch=[0,0,maxploty,  maxploty];
                patch(x_patch, y_patch, 'r','FaceAlpha',.2, 'EdgeColor', 'none');
            end
            xlim([-win(1), win(2)]);
            %get the means / sums, plot on a different figure.
            fr_mean=mean(gau_sdf_store, 1);
            
            %plot
            figure(laser_fr_fig)
            axes('Parent',unit_tab);
            a=subplot(sq2_, sq2_, stype);
            subplot(a)  
            plot(t_plot,fr_mean);
            hold on
            xlabel('Time (s)')
            ylabel('Mean FR (hz) ')
            xlim([-win(1), win(2)]);
            title(laser_labels{stype})      
            maxploty=max(fr_mean);  
            minploty=min(fr_mean);
            if minploty<1
                minploty=0;
            end
            if maxploty<1
                maxploty=2;
            end
            ylim([minploty-1, maxploty+1]);    
            markers=stimMarkers{stype};
            plot( repmat(markers, 1,7), linspace(0, maxploty, 7),'Color', [1 0 0 0.25], 'LineWidth', 1)  %use an odd number in the plotting as laser stamps come in pairs - that way it'll never be square
            set(gca,'FontSize', 14);
            for hh=1:length(markers)/2 %give the laser regions a little shading
                x_patch= [markers((2*hh)-1), markers(2*hh), markers(2*hh), markers(2*hh-1)];
                y_patch=[0,0,maxploty,  maxploty];
                patch(x_patch, y_patch, 'r','FaceAlpha',.2,'EdgeColor', 'none');
            end
            
    end
    
    
end 


%% Smooth firing rate repsonses on the probe by depth, with TTLs. 

laserinf=spikeStruct.laserinf;
laserTTL=laserinf(:,1); %these are the times of laser events
fs=spikeStruct.sample_rate;  %sampling rate
c_channel=spikeStruct.c_channel;  %centre channel for the cluster
min_t=spikeStruct.timeRange(1);  %min, max spikes in the recording
max_t=spikeStruct.timeRange(2);
bl_start=spikeStruct.baseline_st;  %baseline times, if needed.
bl_end=spikeStruct.baseline_end;

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


%% FR plot - one for each cell, plotted in correct depth order.

cellList=1:10; %specifiy which cells to plot (EDIT)
timePlot=[0,1000]; % time to plot, in seconds.
binsize=5; %in seconds

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

TTL_inc=TTLts(TTLts>timePlot(1)& TTLts<timePlot(2));

fr_fig= figure('color','w','NumberTitle','off', 'name',' Unit firing', 'units', 'centimeters', 'pos',[5 3 24 15], 'Color', 'white');
p6 = uipanel('Parent',fr_fig,'BorderType','none'); 

figure(fr_fig)
ADC_plot=subplot('Position',[0.065 0.05 0.88 0.1],  'Parent', p6)

FRsubplots={};
sp_height=0.65/length(cellList);
if sp_height > 0.3
    sp_height=0.24
end


for u=1:length(cellList)
    
    FRsubplots{u}=subplot('Position',[0.065 0.28+((plot_pos2(u)-1)*1.1*sp_height) 0.88 sp_height],  'Parent', p6);
    
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
    spk_count = spk_count(1:end-1);
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
       ylabel(ticklabs{plot_pos2(u)}, 'FontSize', 8 )
    hold on
    
    if plot_pos2(u) >1
        set(gca,'Xticklabel',[]) 
        xlabel('')
    end
    
    aa=gca;
    ymax=aa.YLim(2)
    ymin=aa.YLim(1)
    
    if TTL_inc
       plot(repmat(TTL_inc,1,3), [ymin, ymax, ymax], '-r', 'Color', [1 0 0 0.3], 'LineWidth', 1);
    end
   
end

%now plot bladder pressure and sphincter EMG
subplot(ADC_plot)
hold on
title('Optional ADC plot', 'FontWeight', 'normal')
 xlim([timePlot(1), timePlot(2)])
%now overlay laser TTL markers, if included in the period specified.
if TTL_inc
    ymax=1
    ymin=0
    plot(repmat(TTL_inc,1,3), [ymin, ymax, ymax], '-r', 'Color', [1 0 0 0.2], 'LineWidth', 2);
end
xlim([timePlot(1), timePlot(2)])
title('Optional ADC plot', 'FontWeight', 'normal')
