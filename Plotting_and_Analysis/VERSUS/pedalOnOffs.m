function [pedal_ts] = pedalOnOffs(fn)
%Takes a continuous openEphys ADC file (user must supply the function with the path
%to the file), extracts event times using the rising and falling edges of the signal. 
%NB the threshold for an event (line 25) to be detected might need to be changed depending on 
%the device used to generate the signal.
%Returns a matrix where each row is the start/end of an individual event.

% Read in data
[raw_sig, ADC_ts,  ADC_info] = load_open_ephys_data(fn);  %returns the data, the timebase, and a struct with more useful information

% Set timestamps so that zero is at the start of the recording (as for the
% spiking data from kilosort/PHY)
ADC_tbase=ADC_ts-ADC_ts(1); % timebase should start at zero
fs=ADC_info.header.sampleRate;
%% Detect large gaps in the rec (will generate a warning if the rec is not all in one block
blocks=block_detector(ADC_tbase);

%% Filter (cleans up a bit)

filtlim=20;   %the frequency limit for the LP filter, in Hz (100):
[Db,Da]=butter(2,filtlim/(0.5*fs), 'low' );
ped_sig=filtfilt(Db,Da,raw_sig);  %Our final processed signal.
%% Plot data (shows how events are extracted)

pedThres=3; %threshold for extracting foot down events

figure('Color', 'w', 'units', 'normalized', 'Position', [0.1, 0.2, 0.8, 0.3])
plot(ADC_ts, ped_sig)

%extract times above threshold, detect start and end of periods above
%threshold.
aboveThres=find(ped_sig>pedThres);
sepEvents=find(diff(ADC_ts(aboveThres))>0.5); %find start of separate ped events
starts_=[ADC_ts(aboveThres(1));ADC_ts(aboveThres(sepEvents+1))];
ends_=[ADC_ts(aboveThres(sepEvents));ADC_ts(aboveThres(end))];
pedal_ts=[starts_,ends_];

%plot
for g=1:length(starts_)
    this_event=pedal_ts(g,:);
    p=patch([this_event(1), this_event(2), this_event(2), this_event(1)], [pedThres-0.3,pedThres-0.3, pedThres, pedThres], 'g');
    p.FaceAlpha=0.5;
    p.EdgeColor='none'
end

ylim([-0.5, 3.5])
title('Pedal signal, with individual events marked', 'FontWeight', 'normal')
xlabel('Time (s)')
ylabel('\muV')
end

