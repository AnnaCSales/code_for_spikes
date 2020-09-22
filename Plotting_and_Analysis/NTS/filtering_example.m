
% This is a demo script for analysing mouse ECG data. It does the
% following:
%  - reads in openEphys data for mouse ECG, 
%  - filters it to highlight the parts of the signal relating to the respiratory cycle 
%  - pulls out the number of cycles per second, and calculates the respiratory rate
% Anna Sales, UoB, August 2020.

%To do - write functions for finding resp rate and heart rate
% Go through the events, pull out rate before, during, after etc 
% plot activity around each event.
%% Start by reading in the data and the TTLs

ADC_fn='D:\NTS_Pabitra\NTS_new_data\EMG_ECG\100_ADC2.continuous'; %this is a string containing the location of the data file
[raw_sig, ADC_ts_raw,  ADC_info] = load_open_ephys_data(ADC_fn);  %returns the data, the timebase, and a struct with more useful information

% Read in the events file, for the TTL
TTL_fn='D:\NTS_Pabitra\NTS_new_data\EMG_ECG\all_channels.events'
[events, ev_ts_raw, ev_info] = load_open_ephys_data(TTL_fn);   

%%
analysis_window=[230, 240]; %specify the period we want to analyse, in seconds

% Read in the ECG data:
ADC_tbase=ADC_ts_raw-ADC_ts_raw(1); % timebase should start at zero 
fs=ADC_info.header.sampleRate;  %read the sampling rate


ev_ts=ev_ts_raw-ADC_ts_raw(1); %same time base as for probe data
laser_BNC_chan=4; %where the laser TTL BNC was plugged in to the in/out board
evind=find(events==4);
laser_ts_all=ev_ts(evind);
% ev_ts(2:2:end)=[];

% now cut the data, so that we only work with data in the time range
% specified in 'analysis window'

% the indices of data in the specified range:
keep_=find(ADC_tbase>analysis_window(1) & ADC_tbase<analysis_window(2));

%the timebase and raw signal in the time window:
tbase_win=ADC_tbase(keep_); 
sig_win=raw_sig(keep_);
laser_inc=find(laser_ts_all>analysis_window(1) & laser_ts_all<analysis_window(2));
laser_ts=laser_ts_all(laser_inc);


%% Now do some filtering / editing

zsig=zscore(sig_win); %find the zscored version of the raw signal
mn_sg=mean(zsig); %find the mean of the signal 
std_sig=std(zsig); %find the standard deviation of the signal 

%set anything above 2*stds back to the mean (remove high ECG readings):
nstd=2;
high_ind=find(abs(zsig)>nstd*std_sig);
zsig(high_ind)=0;

%take the moving mean. This will get rid of the higher frequency ECG
%components:
resps=movmean(abs(zsig),8000);

%now high pass filter to remove slow drifts in the baseline:
HPlim=0.2;   %the frequency limit for the HP filter, in Hz:
[Db,Da]=butter(2,HPlim/(0.5*fs), 'high' );
resp_sig=filtfilt(Db,Da,resps);  %Our final processed signal.

%% Now try to filter the ECG signal

% LP filter this time to remove laser artefact / clean up for peak
% detection

HPlim=100;   %the frequency limit for the LP filter, in Hz (100):
[Db,Da]=butter(2,HPlim/(0.5*fs), 'low' );
ecg_sig=filtfilt(Db,Da,sig_win);  %Our final processed signal.

%% Now plot
filt_fig=figure('Color', 'w')
subplot(2,1,1)
plot(tbase_win,resp_sig, 'b');
ylim([-0.2,0.3]);
ylabel('Filtered signal, \muV');
xlim([analysis_window(1), analysis_window(2)]);
xlabel('Time (s)');
hold on
plot([laser_ts, laser_ts], [-0.2, 0.3], 'Color', [1 0 0 0.1])
title('Respiratory signal')


subplot(2,1,2)
plot(tbase_win,ecg_sig, 'b');
ylim([-0.2,0.2]);
ylabel('Filtered signal, \muV');
title('ECG signal')
hold on
xlim([analysis_window(1), analysis_window(2)]);
%%
% find peaks in the data and mark on the plot (may need to play with the
% parameters in the function call)
[pks_resp,times_resp] = findpeaks(resp_sig,fs, 'MinPeakProminence',0.1,'MinPeakDistance',0.1);
times_resp=times_resp+analysis_window(1); %set times relative to start of window
figure(filt_fig)
subplot(2,1,1)
plot(times_resp, pks_resp, 'kX')

% find the negative going peaks in the ECG signal
[pks_ecg,times_ecg] = findpeaks(-1*ecg_sig,fs, 'MinPeakProminence',0.1,'MinPeakDistance',0.1);
times_ecg=times_ecg+analysis_window(1); %set times relative to start of window

subplot(2,1,2)
plot(times_ecg, -1*pks_ecg, 'kX')

%%
 
% work out how many peaks there are, per second:
resp_rate=length(times_resp) ./ range(analysis_window);
 
%print out the final answer:
fprintf('/n Respiratory rate in specified period: %.2f per second /n', resp_rate)
