%plot the 4 bottom channels, to try to see if we did actually record any
%fibres.

%filter settings.
data=[];
[data_raw1, ts_raw, info] = load_open_ephys_data('113_CH17.continuous'); 

fs=info.header.sampleRate; % extract sampling rate
filts=[300,600];  %filter below 40Hz
[Db,Da]=butter(2,filts/(0.5*fs), 'Bandpass' );
  
data(1,:)=filtfilt(Db,Da,data_raw1);  %filtered delta
ts=ts_raw-ts_raw(1);

[data_raw2, ~, ~] = load_open_ephys_data('113_CH1.continuous');
data(2,:)=filtfilt(Db,Da,data_raw2);  %filtered delta

[data_raw3, ~, ~] = load_open_ephys_data('113_CH32.continuous');   
data(3,:)=filtfilt(Db,Da,data_raw3);  %filtered delta

[data_raw4, ~, ~] = load_open_ephys_data('113_CH16.continuous');
data(4,:)=filtfilt(Db,Da,data_raw4);  %filtered delta

[events, ev_ts_raw, ev_info] = load_open_ephys_data('all_channels.events');   
ev_ts=ev_ts_raw-ts_raw(1); %same time base as for probe data
evind=find(events==2);
ev_ts=ev_ts(evind)
ev_ts(2:2:end)=[];
%%
tplot=[0,60];
inds_keep=find(ts>tplot(1) & ts<tplot(2));
ev_inds_keep=find( ev_ts>tplot(1) & ev_ts<tplot(2)) ;
ts_keep=ts(inds_keep);
ylimit=20;
figure('Color', 'w', 'Units', 'Normalized', 'Position', [0.1 0.1 0.8 0.8])
pl=1
for h=4:-1:1
    subplot(4,1,pl);
    plot(ts(inds_keep), data(h,inds_keep));
    xlabel('Time (s)');
    ylabel('\muV');
    hold on
    plot([ev_ts(ev_inds_keep),ev_ts(ev_inds_keep)], [-ylimit ylimit], '--r', 'Color', [1 0 0 ])
    ylim([-ylimit ylimit]);
    if h==1
        title(['Channel:' num2str(h) ' (Lowest)'])
    elseif h==4
        title(['Channel:' num2str(h) ' (Highest)'])
    end
    pl=pl+1;
      
end
%% Pull out spikes on bottom channel
this_data=data(1,inds_keep);

figure('Color', 'w', 'Units', 'Normalized', 'Position', [0.1 0.1 0.8 0.4])
plot(ts_keep, this_data);
xlabel('Time (s)');
ylabel('\muV');
hold on
% plot([ev_ts(ev_inds_keep),ev_ts(ev_inds_keep)], [-ylimit ylimit], 'r')
ylim([-ylimit ylimit]);
xlim(tplot);

% mn_data=nanmean(data(1,:));
% std_data=nanstd(data(1,:));
% nstds=2;  %threshold for cutting - doesn't work because the artefacts are huge
% thres=mn_data+(nstds*std_data);
thres=10;
plot(tplot, [thres, thres], 'k:')
plot(tplot, [-thres, -thres], 'k:')

abv_thr=find(abs(this_data)>thres);
[~,high_inds]=findpeaks(this_data(abv_thr)); %above threshold, at least 2ms apart
abs_high=abv_thr(high_inds); %indicies relative to whole dataset.
% plot(ts_keep(abs_high), thres, 'rs'); %plot to check

% now kill off anything within 5ms of an event, as this is likely to be an
% artefact
ev_times=ev_ts(ev_inds_keep);
% plot(ev_times, thres+1, '*');
to_kill=[];
for t=1:length(ev_times)
    art_inds=find(abs(ts_keep(abs_high)-ev_times(t))<=0.015);
    to_kill=[to_kill, art_inds'];
end

abs_high(to_kill)=[];
% plot(ts_keep(abs_high), thres, 'gs'); %plot to check

%now kill off aretefacts (abnormally high values)

artefacts=find(abs(this_data(abs_high))>13);
% plot( ts(abs_high(artefacts)), 15, 'go');
abs_high(artefacts)=[];
spk_inds=abs_high;

plot(ts(spk_inds), this_data(spk_inds), 'o', 'MarkerSize', 4);

%extract actual spikes.
samp_inc=[40,41]; %number of samples either side of central point.
spk_wins=[spk_inds'-samp_inc(1), spk_inds'+samp_inc(2)];
spikes=[];
for sp=1:size(spk_wins,1)
    spikes(sp,:)=this_data(spk_wins(sp,1):spk_wins(sp,2));
end


figure('Color', 'w')
title('All spikes', 'FontWeight', 'normal')
t_spk=(1:size(spikes,2))/fs;
plot(t_spk*1000, spikes','b')
hold on
plot([t_spk(41),t_spk(41)]*1000, [-15 15], ':k')
xlabel('Time (ms)')
ylabel('\muV')
title('All extracted spikes', 'FontWeight', 'normal')