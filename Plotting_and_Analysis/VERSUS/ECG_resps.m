% This is a script for analysing rat ECG/respiratory data from the VERSUS prep. It does the
% following:
%  - reads in openEphys data for  ECG,
%  - filters it to highlight the parts of the signal relating to the respiratory cycle / ECG

% You need to update the filename fields, e.g. ADC_fn and TTL_fn, so that
% the script can find the ADC data and the TTL events, and also make sure
% that the 'laser_BNC_chan' correctly shows where the laser BNC was plugged
% in to the input/output board.
% Anna Sales, UoB, August 2020.


%% Start by reading in the raw ADC data and the TTL information


ADC_fn=[pwd '\100_34_2.continuous']; %this is a string containing the location of the data file
[raw_sig, ADC_ts_raw,  ADC_info] = load_open_ephys_data(ADC_fn);  %returns the data, the timebase, and a struct with more useful information

% Set timestamps so that zero is at the start of the recording (as for the
% spiking data from kilosort/PHY)
ADC_tbase=ADC_ts_raw-ADC_ts_raw(1); % timebase should start at zero
fs=ADC_info.header.sampleRate;  %read the sampling rate

%% Detect large gaps in the rec

blocks=block_detector(ADC_tbase);

%% optional - specify a time period (otherwise use entire recording - comment out as req)

% analysis_window=[230, 240]; %specify the period we want to analyse, in seconds

% now cut the data, so that we only work with data in the time range
% specified in 'analysis window'

% the indices of data in the specified range:


%the timebase and raw signal in the time window:
ADC_ts=ADC_tbase;
raw_sig=raw_sig(keep_);

%% Filter, find peaks representing heart beat/ breathing

[ecg, resps]=return_filtered_data(raw_sig, fs);

% find peaks in the data and mark on the plot (may need to play with the
% parameters in the function call)
[pks_resp,times_resp] = findpeaks(resps,fs, 'MinPeakProminence',0.08,'MinPeakDistance',0.1);
times_resp=times_resp+analysis_window(1); %set times relative to start of window

% store a continuous log of gaps between breaths
resp_gaps=[0; diff(times_resp)];
times_gaps=[times_resp-[0;0.5*diff(times_resp)]];  %time points for these gaps (halfway between breaths)

% find the negative going peaks in the ECG signal
[pks_ecg,times_ecg] = findpeaks(-1*ecg,fs, 'MinPeakProminence',200,'MinPeakDistance',0.05);
times_ecg=times_ecg+analysis_window(1); %set times relative to start of window

% store the continuous rate (per min) as 60/(time between consecutive beats)
ecg_rate=[0;60./diff(times_ecg)];

%% Now plot 20 secs of data (helps to sanity check, and
% fine tune parameters for peak finding above)

filt_fig=figure('Color', 'w')

subplot(3,1,1)
in_win=find(ADC_ts>analysis_window(1) & ADC_ts<analysis_window(1)+10);
plot(ADC_ts(in_win),raw_sig(in_win), 'k');
ylim([-500,300]);
ylabel('Raw signal, \muV');
xlim([analysis_window(1), analysis_window(1)+10]);
xlabel('Time (s)');
hold on
plot(times_resp, pks_resp, 'kX')
title('Example data: raw signal', 'FontWeight', 'normal')


subplot(3,1,2)
in_win=find(ADC_ts>analysis_window(1) & ADC_ts<analysis_window(1)+10);
plot(ADC_ts(in_win),resps(in_win), 'b');
ylim([-0.1,0.15]);
ylabel('Filtered signal, \muV');
xlim([analysis_window(1), analysis_window(1)+10]);
xlabel('Time (s)');
hold on
plot(times_resp, pks_resp, 'kX')
title('Example data: respiratory signal', 'FontWeight', 'normal')


subplot(3,1,3)
plot(ADC_ts(in_win),ecg(in_win), 'b');
ylim([-500,500]);
ylabel('Filtered signal, \muV');
title('ECG signal','FontWeight', 'normal')
hold on
plot(times_ecg, -1*pks_ecg, 'kX')
xlim([analysis_window(1), analysis_window(1)+10]);

%% Plot firing rate vs ECG

cellList=[2]
timePlot=[10, 30]; % time to plot, in seconds.
binsize=0.1; %in seconds
fs=spikeStruct.sample_rate;  %sampling rate
c_channel=spikeStruct.c_channel;  %centre channel for the cluster
ts_= spikeStruct.timesSorted{iUnit}; 
ts_ind=find(ts_>=timePlot(1)& ts_<=timePlot(2));
ts_window=ts_(ts_ind);  %store all the data that's been cut.
isr=1./diff(ts_window);
tbin_edges = timePlot(1):binsize:timePlot(2);
tbin_centers = tbin_edges(1:end-1)+binsize/2;
t_plot=tbin_centers;


%bin and histogram
spk_count = histc(ts_window,tbin_edges);
spk_count = spk_count(1:end-1) ./binsize;
%plot
figure('Color', 'w')
subplot(2,1,1)
% plot(ts_window(2:end), isr)
mybar=bar(t_plot, spk_count);
mybar.FaceColor=[0 0 1];
mybar.FaceAlpha=0.7;
mybar.LineWidth=1;
xlabel('Time (s)');
xlim(timePlot);
ylabel('FR (Hz)')

subplot(2,1,2)
in_win=find(ADC_ts>timePlot(1) & ADC_ts<timePlot(2));
plot(ADC_ts(in_win),ones(length(in_win), 1), '.b');
% ylim([-500,500]);
ylabel('Filtered signal, \muV');
title('ECG signal','FontWeight', 'normal')
hold on
% plot(times_ecg, -1*pks_ecg, 'kX')
% xlim([analysis_window(1), analysis_window(1)+10]);

