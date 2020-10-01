% This script contains a variety of useful commands for doing an initial look at data from openEphys
%Comment \ uncomment \ run by highlighting whatever you need. 
% Anna Sales April 2017
%% 
[event_data_all, event_ts_all, event_info] = load_open_ephys_data('D:\Versus\290920\290920\Rec5\all_channels_5.events');
[data_all, ts_offset, info] = load_open_ephys_data('D:\Versus\290920\290920\Rec5\126_CH17_5.continuous');   
%[spike_data, spike_ts, spike_info] = load_open_ephys_data('C:\openEphys\Recordings\2017-02-03_12-18-05\SE0.spikes');
fs=info.header.sampleRate; % extract sampling rate

ts_all=ts_offset-ts_offset(1);

win=[0,1200];
keep_inds=find(ts_all>win(1) & ts_all<win(2));

data=data_all(keep_inds);

event_keep_inds=find(event_ts>win(1) & event_ts<win(2));
event_data=event_data_all(event_keep_inds);
event_ts=event_ts_all(event_keep_inds)-ts_offset(1);



fs=info.header.sampleRate; % extract sampling rate
figure
title('raw data')
plot(ts, data) 
xlabel('Time(s)');
ylabel('Voltage (\mu V)');
% ylim([-400, 400]);
% xlim([ts2(1), ts2(end)]);
hold on
plot(event_ts(event_data==2), 100, 'ro');
plot(event_ts(event_data==6), 100, 'gs');
%block_times=cumsum([110.0, 84+244.0, 108+571.0, 546.0+55, 113+280.0, 348.0+104, 128.0] );
%hold on;
%plot(repmat(block_times, 50, 1), linspace(-1000, 1000, 50) , '.')
%%

[data_all1, ts_offset, ~] = load_open_ephys_data('D:\Versus\290920\290920\Rec5\126_CH25_5.continuous', 'Indices',keep_inds(1):keep_inds(end) );   
[data_all2, ~, ~] = load_open_ephys_data('D:\Versus\290920\290920\Rec5\126_CH14_5.continuous', 'Indices',keep_inds(1):keep_inds(end) );   
[data_all3, ~, ~] = load_open_ephys_data('D:\Versus\290920\290920\Rec5\126_CH30_5.continuous', 'Indices',keep_inds(1):keep_inds(end) );   
[data_all4, ~, info] = load_open_ephys_data('D:\Versus\290920\290920\Rec5\126_CH3_5.continuous', 'Indices',keep_inds(1):keep_inds(end) );   

%%
%Filter data - COMMENT OUT IF NOT NEEDED

filterRange=[300 6000];    %freq range to keep
[b,a]=butter(3,[filterRange(1)/(0.5*fs) filterRange(2)/(0.5*fs)], 'bandpass');  %Butterworth filter, with freq limits converted to pi-rads/sample
data=filtfilt(b,a,data);  %filtered data
% figure
% plot(ts, data) 
% hold on;
% xlabel('Time(s)');
% ylabel('Voltage (\mu V)');
% title('Filtered data')



%% 
%Get rid of extreme values
cutoff=1000;

cut_inds=find(data>cutoff | data<(-1*cutoff));  
data(cut_inds)=[];
ts(cut_inds)=[];
figure
plot(ts/60, data) 
hold on;
xlabel('Time(m)');
ylabel('Voltage (\mu V)');

%% 

%plot the timestamps of spike events - to see where the data is and where
%the periods of silence are. 
% figure;
% plot(spike_ts/60,ones(length(spike_ts),1),'.') % Raster (minutes)   %timestamps for each spike - helps to see where recording starts and ends.
% hold on
% 


%% 
 
 %find the timestamps of laser event (which have event data '2' - check which number corresponds to which event on your data)
 % eventData 0 to 7 corresponds to TTLs 1 to 8
 %NB '0' corresponds to TTL input 1 OR to network events. If you have have
 %used TTL 1 will have to distinguish between the two types using
 %eventInfo.eventType. In this file '5' denotes a network event whilst '3'
 %denotes a TTL.#

laser_ts=event_ts( find(event_data==6));  %find of elements which have event = 2
message_ts=event_ts(find(event_data==0)); %find message timestamps
%% 
%Section to look only at one period 
t_start=3000    ;  %NB ts only recorded during actual periods - nothing in between.
t_end=4000;
    
[~, start_ind]=min(abs(ts-t_start));
[~, end_ind]=min(abs(ts-t_end));

[~, start_laser]=min(abs(laser_ts-t_start));
[~, end_laser]=min(abs(laser_ts-t_end))  ;
laser_ts2=laser_ts(start_laser:end_laser);
laser_ts2(laser_ts2>t_end)=[];
laser_ts2(laser_ts2<t_start)=[];

data2=data(start_ind:end_ind); 
ts2=ts(start_ind:end_ind);

num=size(data2, 1);

laser_marks=linspace(-250, 200, 50); %for plotting laser line;


figure

    hold on
     plot(ts2, data2);
    xlabel('Time(s)');
    ylabel('Voltage (\mu V)');
    ylim([-400, 400]);
    xlim([ts2(1), ts2(end)]);
    for p=1:size(laser_ts2)
        l=plot(laser_ts2(p)*ones(50, 1), laser_marks(1,:), 'r');
        l.Color(4)=0.8;
        l.LineWidth=1;
    end
 
  
fp=[94712100  96678000 98551500 191753700 234071400 235430100 237419700  237889500 ]  ;
% fp=[93268500, 98431800, 106696924, 127950000, 130680000, 133080000];

%for 170517 - 4th one is clon
% fp=[94712100,  96678000, 98551500, 191753700, 234071400, 235430100 ,237419700 , 237889500 ]

fp=fp./(fs);  %put into seconds. 

   
    fp_marks=linspace(-300, 300, 50); %for plotting fp line;
    for p=1:size(fp, 2)
        fp(p);
        ff=plot(fp(p).*ones(50, 1), fp_marks, 'g', 'LineWidth', 1);
ff.Color(4)=0.5;
    end
       
    
%%
%Now specify a region of interest in seconds and display raw data on multiple plots with
%laser timestamps. 100 s per plot.

t1=2950;
n_plots=1;
win=30; %time window in s
close all;
for tt=1:n_plots
  
    t_start=t1+((tt-1)*win);  %NB ts only recorded during actual periods - nothing in between.
    t_end=t_start+win;

    [~, start_ind]=min(abs(ts-t_start));
    [~, end_ind]=min(abs(ts-t_end));
    
    [~, start_laser]=min(abs(laser_ts-t_start));
    [~, end_laser]=min(abs(laser_ts-t_end))  ;

    laser_ts2=laser_ts(start_laser:end_laser);
    data2=data(start_ind:end_ind); 
    ts2=ts(start_ind:end_ind);

    num=size(data2, 1);
    laser_marks=linspace(-250, 200, 50); %for plotting laser line;


    figure  
   plot(ts2, data2);
    hold on;
    xlabel('Time(s)');
    ylabel('Voltage (\mu V)');
    ylim([-400, 400]);
    xlim([ts2(1), ts2(end)]);
    
    for p=1:size(laser_ts2)
        l=plot(laser_ts2(p)*ones(50, 1), laser_marks(1,:), 'r');
        l.Color(4)=0.7;
        
    end
   
    
    %add in a manually generated vector of fp times, if needed - save a vector
    %fp into memory. Extract this info from the messages file.
    %for 7/2/17 recordings:

%     fp=[233442000,237317700,237390300, 334780500, 342569700,356349300,356759400,367151400,367524900];
%     fp=fp./(fs);  %put into minutes. 

   
%     fp_marks=linspace(-300, 300, 50); %for plotting laser line;
%     for p=1:size(fp, 2)
%         fp(p);
%         plot(fp(p).*ones(50, 1), fp_marks, 'g');
% 
%     end
end

