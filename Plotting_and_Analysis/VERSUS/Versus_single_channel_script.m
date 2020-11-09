%Quick script to read in openEphys continuous data files (just the raw
%data), filter it, cut out extreme values, extract spikes. Produces a
%vector of spike times and a matrix of spikes, windowed around the peak (40
%samples). Follows Quiroga (2004, Neural Computation) for details of
%threshold and definitions of noise. 


% Allows basic spike sorting based on the data from one channel only, using
% spike amplitude, width, and projections along principal components.

%Anna Sales October 2020
% WHERE IS YOUR DATA?
data_dir='C:\Versus data\281020\2020-10-28_16-13-14\CAR\';

% raw data
[data_raw, ts_raw, info] = load_open_ephys_data([data_dir '102_CH19.continuous']);   
first_timepoint=ts_raw(1);
ts_raw=ts_raw-first_timepoint;

[events, ts_ev, info_ev] = load_open_ephys_data([data_dir 'all_channels.events']);   
ts_ev=ts_ev-first_timepoint;

% events
% TTLs=returnTTLs(data_dir);  %returns manual and digital TTLs

%footshocks
fsts=footshock_TTLs;
% fsts_=TTLs.digital{7};
% fsts_(2:2:end)=[];  %they come in pairs. Just mark the onset.


%pedal
% pedal_ts=spikeStruct.TTLs.digital{7};
% sep_events_inds=find(diff(pedal_ts)>0.2); %separate pedal pushes
% pedts_=[pedal_ts(1);pedal_ts(sep_events_inds+1)];


%%  select a time window of interest (optional - just use entire rec if needed)

t_start=0;  %time window bounds, in seconds. 
t_end=4000;

t_range=[t_start, t_end];
t_keep=find(ts_raw>t_start & ts_raw<=t_end); %indices within selected time range

data=data_raw(t_keep);   %cut data and timestamps according to time selected.
ts=ts_raw(t_keep);
fs=info.header.sampleRate; % extract sampling rate

% cut events to match the same range
fsts=fsts_(find(fsts_>t_start & fsts_<=t_end));
% pedts=pedts_(find(pedts_>t_start & pedts_<=t_end));
%% Filter data (if not already done)

% filterRange=[300 6000];    %freq range to keep
% [b,a]=butter(3,[filterRange(1)/(0.5*fs) filterRange(2)/(0.5*fs)], 'bandpass');  %Butterworth filter, with freq limits converted to pi-rads/sample
% data_=filtfilt(b,a,data_);  %filtered data

%%  Get rid of the worst extreme values
cutoff=1000;
cut_inds=find(data>cutoff | data<(-1*cutoff));  
data(cut_inds)=[];
ts(cut_inds)=[];

%% Now estimate the noise and set a threshold for spike detection

times_noise=5;  
%The number of  from the mean that we're using as a
%threshold for spike extraction

sd_est=median(abs(data)./0.6745) ;   %estimate of the noise, (following Quiroga(2004) )
% or just use standard deviations
% sd_est= std(abs(data));
threshold=-times_noise*sd_est  %negative threshold

%% Display the filtered data and the thresholds for spike extraction
figure('Color', 'w')
plot(ts, data, 'Color', [0 0 1 0.4])
hold on;
xlabel('Time(s)');
ylabel('Voltage (\mu V)');
% plot(t_range, [threshold, threshold], 'g:', 'LineWidth', 1.5)
plot(t_range, [threshold, threshold], 'g:', 'LineWidth', 1.5)
ylim([-200, 200])
title('Filtered data with spike threshold shown', 'FontWeight', 'normal')
%%  Extract spikes.
%Look only at points above the threshold, find peaks  in the data, extract
%a window around each threshold crossing.

inds_above=find(data<threshold);
data_above=data(inds_above);
[pks,locs] = findpeaks(-1.*data_above);  %find locations of minima

%now pull out the actual spikes.
peak_inds=inds_above(locs);  %index of spikes in full dataset
num_spikes=length(peak_inds);

%Pull out a window of data for each spike, centred on the peak and 40 samples
%long. As in openEphys spikesorter, centre the peak 0.3ms into the window

spike_data=zeros(num_spikes, 80);

for j=1:num_spikes
    spike_data(j,:)=data(peak_inds(j)-39:peak_inds(j)+40);
end

spike_times=ts(peak_inds);  %spike timestamps.
spike_tbase=((1:80)./fs).*1000; %spike timebase in ms
%%  Look at the spikes, before clean up
spkfig=figure('Color', 'w','Units', 'normalized', 'Position', [0.1 0.1 0.62 0.4])
subplot(1,2,1)
plot(spike_tbase, spike_data)
hold on
title('Plot of all spikes extracted, before clean-up', 'FontWeight', 'normal')
xlim([0, spike_tbase(end)])
xlabel('ms')
ylabel('\mu V')

%% Clean up floating voltages (and wrongly windowed spikes), before and after the spike centre
[r_noise, ~]=find(abs(spike_data)>200);

sp_early=spike_data(:, 1:32);
[r_e,~]=find(abs(sp_early)>40);  %base the limit on 8*estimate of noise, above.

sp_late=spike_data(:, 48:end);
[r_l,~]=find(abs(sp_late)>40);

sp_med=spike_data(:, 39:41);
[r_m,~]=find((sp_med)>0) ;  %points at this point should be negative

%artefacts - amend as appropriate for the recording **TAKE THIS BACK
%OUT!**
% sp_val=spike_data(:, 8:9);
% [r_valley, ~]=find((sp_val)<-300);
% killrows=unique([r_e', r_l', r_m', r_valley']); %rows to delete

killrows=unique([r_e', r_l', r_m', r_noise']); %rows to delete

spike_times(killrows, :)=[];
spike_data(killrows, :)=[];

%% 
spikes=struct;
spikes.times=spike_times;
spikes.data=spike_data;

figure(spkfig)
subplot(1,2,2)
plot(spike_tbase, spike_data)
hold on
title('Plot of spikes, cleaned', 'FontWeight', 'normal')
xlabel('ms')
ylabel('\mu V')
xlim([0, spike_tbase(end)])

% figure('Color', 'w', 'units', 'normalized', 'Position', [0.1 0.1 0.6, 0.3])
% subplot(2,1,1)
% plot([spike_times, spike_times], [-1,1], 'k')
% ylim([-3,3])
% xlim([t_range(1), t_range(2)])
% xlabel('Time (s)')
% title('All spikes, raster plot', 'FontWeight', 'normal','FontSize', 10)
% set(gca, 'FontSize', 8)
% hold on
% plot([fsts, fsts], [-3,3], 'Color', [1 0 0 0.2])
% % plot([pedts, pedts], [-3,3], 'Color', [0 1 0.1 0.2])

% subplot(2,1,2)
% binwin=1;
% binedges=t_range(1):1:t_range(2);
% bincenters=binedges-binwin/2;
% bincenters(1)=[];
% 
% spk_rates=histcounts(spike_times, binedges);
% plot(bincenters, spk_rates)
% xlim([t_range(1), t_range(2)])
% xlabel('Time (s)')
% ylabel('Spk rate (Hz)')
%% Now extract tables of features - amplitude, width, projections along top principal components.
% extract peak amplitude and FWHM.
plotExamples=1;  % will plot 5 examples if =1
[amps,fwhms]=find_wf_features_matrix(spike_data, fs, plotExamples); 

% kill off any for which it was not possible to define fwhm etc (should be
% zero, but just in case)

notDefined=find(isnan(fwhms));
if notDefined
    spike_data(notDefined,:)=[];
    amps(notDefined)=[];
    fwhms(notDefined)=[];
    spike_times(notDefined)=[];
end

[coeff,score,latent] = pca(spike_data);
% projections along first and second PC:
scorefirst=score(:, 1);
scoresecond=score(:,2); 
%%

% Now take a look in feature space
figure('Color', 'w', 'units', 'normalized', 'Position', [0.1, 0.1, 0.6, 0.4])
subplot(1,3,1)
scatter(amps, fwhms, 'ro')
xlabel('Amplitude (\muV)')
ylabel('FWHM (ms)')

subplot(1,3,2)
scatter(scorefirst, scoresecond, 'bo')
xlabel('Projection along PC1')
ylabel('Projection along PC2')

subplot(1,3,3)
scatter(scorefirst,fwhms, 'go')
xlabel('Projection along PC1')
ylabel('FWHM (ms)')

%% Do some kmeans clustering using a projection that works
% featMatrix=[scorefirst,scoresecond];
featMatrix=[amps, fwhms];

%try to determine optimum number of clusts
total_dist=[];
for k=1:4
      [~,~,dists]=kmeans(featMatrix, k);
     total_dist(k)=sum(dists);
end
figure('Color', 'w')
plot(1:4, total_dist)
xlabel('# clusts')
ylabel('Total distance, points from allocated centroid')
%%
numclusts=2;
if numclusts==2
    [idx, c]=kmeans(featMatrix, numclusts);
    
    gp1=find(idx==1);
    gp2=find(idx==2);
    gp3=NaN;
    figure('Color', 'w')
    scatter(featMatrix(:,1),featMatrix(:,2), 'k.')
%     xlabel('Projection along PC1')
%     ylabel('Projection along PC2')
    xlabel('Amplitude (\muV)')
    ylabel('FWHM (sm)')

    hold on
    scatter(featMatrix(gp1,1),featMatrix(gp1,2), 'ro')
    scatter(featMatrix(gp2,1),featMatrix(gp2,2), 'bo')
elseif numclusts==3
    [idx, c]=kmeans(featMatrix, numclusts);
    
    gp1=find(idx==1);
    gp2=find(idx==2);
    gp3=find(idx==3);
    figure('Color', 'w')
    scatter(featMatrix(:,1),featMatrix(:,2), 'k.')
    xlabel('Amplitude (\muV)')
    ylabel('FWHM (sm)')
    hold on
    scatter(featMatrix(gp1,1),featMatrix(gp1,2), 'ro')
    scatter(featMatrix(gp2,1),featMatrix(gp2,2), 'bo')
    scatter(featMatrix(gp3,1),featMatrix(gp3,2), 'go')
end
%% Now plot raster of spiking for either unit
spkGroup1=spike_times(gp1);
spkGroup2=spike_times(gp2);

figure('Color', 'w', 'units', 'normalized', 'Position', [0.1 0.1 0.6, 0.8])

subplot(6,1,1)
plot([spkGroup1, spkGroup1], [-1,1], 'k')
ylim([-3,3])
xlim([t_range(1), t_range(2)])
xlabel('Time (s)')
title('Cluster #1, from single channel data', 'FontWeight', 'normal','FontSize', 8)
set(gca, 'FontSize', 8)
hold on
plot([fsts, fsts], [-3,3], 'Color', [1 0 0 0.2])
plot([pedts, pedts], [-3,3], 'Color', [0 1 0.1 0.2])



subplot(6,1,2)
plot([spkGroup2, spkGroup2], [-1,1], 'k')
ylim([-3,3])
xlim([t_range(1), t_range(2)])
title('Cluster #2, from single channel data', 'FontWeight', 'normal','FontSize', 8)
xlabel('Time (s)')
set(gca, 'FontSize', 8)
hold on
plot([fsts, fsts], [-3,3], 'Color', [1 0 0 0.2])
plot([pedts, pedts], [-3,3], 'Color', [0 1 0.1 0.2])


subplot(6,1,3)
if ~isnan(gp3)
    spkGroup3=spike_times(gp3);
    plot([spkGroup3, spkGroup3], [-1,1], 'k')
    ylim([-3,3])
    xlim([t_range(1), t_range(2)])
    title('Cluster #3, from single channel data', 'FontWeight', 'normal','FontSize', 8)
    xlabel('Time (s)')
    set(gca, 'FontSize', 8)
    hold on
    plot([fsts, fsts], [-3,3], 'Color', [1 0 0 0.2])
    plot([pedts, pedts], [-3,3], 'Color', [0 1 0.1 0.2])

else
    spkGroup3=[];
    title('NO GROUP 3 FROM SINGLE CHAN', 'FontWeight', 'normal', 'FontSize', 8)
end


subplot(6,1,4)
spks=spikeStruct.timesSorted{6};
plot([spks, spks], [-1,1], 'k')
ylim([-3,3])
xlim([t_range(1), t_range(2)])
xlabel('Time (s)')
title('Cluster #6, from multiunit data', 'FontWeight', 'normal', 'FontSize', 8)
set(gca, 'FontSize', 8)
hold on
plot([fsts, fsts], [-3,3], 'Color', [1 0 0 0.2])
plot([pedts, pedts], [-3,3], 'Color', [0 1 0.1 0.2])

subplot(6,1,5)
spks=spikeStruct.timesSorted{3};
plot([spks, spks], [-1,1], 'k')
ylim([-3,3])
xlim([t_range(1), t_range(2)])
xlabel('Time (s)')
title('Cluster #3, from multiunit data', 'FontWeight', 'normal','FontSize', 8)
set(gca, 'FontSize', 8)
hold on
plot([fsts, fsts], [-3,3], 'Color', [1 0 0 0.2])
plot([pedts, pedts], [-3,3], 'Color', [0 1 0.1 0.2])


subplot(6,1,6)
spks=spikeStruct.timesSorted{2};
plot([spks, spks], [-1,1], 'k')
ylim([-3,3])
xlim([t_range(1), t_range(2)])
xlabel('Time (s)')
title('Cluster #2, from multiunit data', 'FontWeight', 'normal','FontSize', 8)
set(gca, 'FontSize', 8)
hold on
plot([fsts, fsts], [-3,3], 'Color', [1 0 0 0.2])
plot([pedts, pedts], [-3,3], 'Color', [0 1 0.1 0.2])

%% Now plot the auto and xcors from the single chan data

% Make a fake spikeStruct
spikeStruct_single.timesSorted{1}=spkGroup1;
spikeStruct_single.timesSorted{2}=spkGroup2;
spikeStruct_single.timesSorted{3}=spkGroup3;

spikeStruct_single.timeRange=[0, 1200];
spikeStruct_single.sample_rate=fs;

pmssingle.cell_list=1:2;
pmssingle.binsize=1;  %binsize, in ms
pmssingle.window=300;  %window size, in ms.
pmssingle.baseline = [] ; %leave empty if no baseline!

[corfig] = cell_xcors(spikeStruct_single, pmssingle)

%% Plot xcors, to check isolation - NB don't plot too many as it's slow
pms2.cell_list=[2,3,6];
pms2.binsize=1;  %binsize, in ms
pms2.window=300;  %window size, in ms.
pms2.baseline = [] ; %leave empty if no baseline!

[corfig] = cell_xcors(spikeStruct, pms2)