

[data_raw1, ts_raw1, info1] = load_open_ephys_data('113_CH1.continuous');   
filts=[300,600]   %filter below 40Hz
[Db,Da]=butter(2,filts/(0.5*fs), 'Bandpass' );
data_1=filtfilt(Db,Da,data_raw1);  %filtered delta

[data_raw2, ts_raw2, info2] = load_open_ephys_data('113_CH17.continuous');
data_2=filtfilt(Db,Da,data_raw2);  %filtered delta

[events1, ev_ts, ev_info] = load_open_ephys_data('all_channels.events');   
fsind=find(events1==2);

fs=info1.header.sampleRate; % extract sampling rate
%%
figure('Color', 'w')
subplot(2,1,1)
plot(ts_raw1, data_2)
hold on
plot([ev_ts(fsind),ev_ts(fsind)], [-500 500], 'r')
ylim([-50,50])
subplot(2,1,2)
plot(ts_raw1, data_1)
hold on
plot([ev_ts(fsind),ev_ts(fsind)], [-500 500], 'r')
ylim([-50,50])
%%
LASER_BNC_CH=7;
laser_ts=unique(ev_ts( find(events1==LASER_BNC_CH)));  %find of elements which have laser timestamp.

% FS=1;
% fs_ts=unique(ev_ts( find(events1==FS)));  %find of elements which have laser timestamp.
%% 
%select a time window of interest 
% 

t_start=0; 
t_end=1000;

t_range=[t_start, t_end];

t_keep1=find(ts_raw1>t_start & ts_raw1 <t_end);
data1=data_raw1(t_keep1);   %cut data according to time selected.
ts1=ts_raw1(t_keep1);

t_keep2=find(ts_raw2>t_start & ts_raw2 <t_end);
data2=data_raw2(t_keep2);   %cut data according to time selected.
ts2=ts_raw2(t_keep2);

laser_keep=find(laser_ts>t_start & laser_ts <t_end);
laser_ts1=laser_ts(laser_keep)
%%

yrange=[min(data1), max(data1)];

figure
subplot(2,1,1)
title('Raw data, zeroed, laser timetamps overlaid(r=raw)')
hold on
ylabel('\muV')
xlabel('time (s)')
t_first=ts1(1);
plot(ts1-t_first, data1)
plot(repmat(laser_ts1-t_first, 1,2), [yrange(1), yrange(2)], '--','Color', [1 0 0 0.3])

subplot(2,1,2)
hold on
ylabel('\muV')
xlabel('time (s)')
t_first=ts1(1);
plot(ts2-t_first, data2)
plot(repmat(laser_ts1-t_first, 1,2), [yrange(1), yrange(2)], '--','Color', [1 0 0 0.3])


%%


for g=1:length(laser_pulses)
%     laser_=laser_pulses{g};
laser_=ttt{g};
    xs_=length(laser_);
    plot(laser_, 2*ones(xs_,1), 'cx', 'MarkerSize', 5)
end
    


for g=1:length(multi_stims)
laser_=multi_stims{g}(:,1);
xs_=length(laser_);
    plot(laser_, 10*ones(xs_,1), 'mx', 'MarkerSize', 7)
end

    
    
    
    
    
%     
% subplot(2,1,2)
% plot(ts2, data2)
% hold on
% ylabel('\muV')
% xlabel('time (s)')
% plot(repmat(laser_ts1, 1,2), [yrange(1), yrange(2)], 'Color', [1 0 0 0.3])