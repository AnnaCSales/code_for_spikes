   
[data, ts, info]=load_open_ephys_data('D:\Versus\0408\Recording1\FGPA\100_CH10_3.continuous');
fs=info.header.sampleRate;
ts=ts-ts(1);
%Filter data for SWS (makes the above filtering options pointless)
rangeDelta=[0.5 4];   %filter below 2hz
[Db,Da]=butter(2,[rangeDelta(2)/(0.5*fs) ]);

EEG_delta=filtfilt(Db,Da,data);  %filtered delta
EEG_delta_p=abs(hilbert(EEG_delta));  %power in the delta
delta_phase = (angle(hilbert(EEG_delta)))';  %phase of the delta
%% Assign every time point and TTL into one of 4 phase bins.

[N,edges, binIDs]=histcounts(delta_phase, [-pi, -pi/2, 0, pi/2, pi]);
pb_labels={'-\pi to -\pi/2 (just after trough) ' '-\pi/2 to 0 (approaching peak)' '0 to \pi/2 (just after peak)' '\pi/2 to pi (approaching trough)'}

% [N,edges, binIDs]=histcounts(abs(delta_phase), [0,  pi/2, pi]);
% pb_labels={'Peak region' 'Trough region'}
TTLchan=3;  %this is the channel to take TTLs from
TTLts=spikeStruct.TTLs.digital{TTLchan};
TTLts(2:2:end)=[];

phasebin=[];
for y=1:length(TTLts)
    TTL_=TTLts(y);
    ind_in_data=find(ts==TTL_);
    phasebin(y)=binIDs(ind_in_data);
end

%% Plot example figure showing phase and filtered LFP

nTTLs_inc=3
plotlims=find(ts>TTLts(1) & ts<TTLts(nTTLs_inc));
% plotlims=[1:3e5]; % time period to plot.

figure('Color', 'w')
subplot(3,1,1)
plot(ts(plotlims), EEG_delta(plotlims))
xlabel('Time (s)')
ylabel({'LFP in delta' 'range (0.5-4Hz), \muV'})

subplot(3,1,2)
plot(ts(plotlims), delta_phase(plotlims))
yticks([-2*pi, -pi, 0 pi, 2*pi])
yticklabels({'-2\pi' '-\pi' '0' '\pi', '2pi'})
xlabel('Time (s)')
ylabel('Phase angle (rads)')

% Peak is at zero, trough is at pi, rising edge is -pi back up to zero.
 subplot(3,1,3)
 plot(ts(plotlims), binIDs(plotlims), 'ro')
 xlabel('Time (s)')
 ylabel('Bin assignment')
 hold on
 plot([TTLts(1:nTTLs_inc) TTLts(1:nTTLs_inc)], [1,length(pb_labels)], 'g')

%% Plotting activity around a TTL - change fields as required
% Bin 1 Trough to halfway up
% Bin 2 Just before peak
% Bin 3 Just after peak
% Bin 4 Approaching trough 

win=[0.25, 0.5] ; %specify a window, time before and after event to consider for pinches
binwin=0.01;
foot_sp_fig      = figure('color','w','NumberTitle','off', 'name','Spiking around footshock / pinch TTLs', 'units', 'centimeters', 'pos',[5 2 24 17]);
foot_spTabGroup = uitabgroup(foot_sp_fig,'TabLocation','Left');

for iUnit=1:nclusts

  figure(foot_sp_fig)
  footunit_tab_sp = uitab(foot_spTabGroup, 'Title', ['Cluster ' num2str(iUnit)],'BackgroundColor','w');
  axes('Parent',footunit_tab_sp);

  %pull out the relevant spike times
  ts_= spikeStruct.timesSorted{iUnit};

  event_ts=[];
  spk_count_all=[];
  
  for pb=1:length(pb_labels)
      
      TTL_to_plot=TTLts(phasebin==pb);  %update this as needed. These are the TTLs to plot.
      nTTL=length(TTL_to_plot);
      
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

      end
 
  %plot mean firing rate
  subplot(4,1,pb)
  binned_mean_fr=mean(spk_count_all,1);
  bar(t_plot,binned_mean_fr);
  hold on
  xlabel('Time (s)')
  ylabel('Mean spike count')
  text(0.3, 2,[num2str(1000*binwin) 'ms bins'])
  xlim([-win(1), win(2)]);
  maxploty=max(binned_mean_fr)+1;
  minploty=0;
  plot( zeros(1,7), linspace(0, maxploty, 7), 'r', 'LineWidth', 1)  %use an odd number in the plotting as laser stamps come in pairs - that way it'll never be square
  ylim([0,3])
  set(gca, 'FontSize', 11);
  title(pb_labels{pb}, 'FontWeight', 'normal')
  end
  
end
