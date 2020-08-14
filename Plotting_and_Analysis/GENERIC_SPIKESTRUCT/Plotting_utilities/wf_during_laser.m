function [wf_fig] = wf_during_laser(spikeStruct,laser_events, duration, fn)
%plots waveforms during laser stimulation, and at other times, for
%comparison. 
% Expects - spikeStruct, list of laser times, duration of each pulse (in ms),
% filename of raw data file.
% Returns - handle to figure generated.

wfparams.dataDir=[pwd '\'];; %where is raw data
wfparams.fileName=fn;
wfparams.nCh = 32;
wfparams.dataType='int16';
wfparams.wfWin = [-40 41];  %number of samples before and after spike peak
wfparams.nWf = 2000;   %number of waveforms to return (if they are there)
wfparams.nBad=0; %nobad channels.
fs=spikeStruct.sample_rate;
%make a table of times to consider,e.g during TTLs
spikes_duringTTL=[];
clusterIDs_duringTTL=[];  
all_spk_ts=spikeStruct.st;
all_spk_IDs=spikeStruct.clu;
duration_s=duration/1000;

for p=1:length(laser_events)   %get a list of spikes / cluster IDs which happened during laser.
    event_time=laser_events(p);
    inds_duringTTL=find(all_spk_ts>event_time & all_spk_ts< (event_time+duration_s) ); 
    spikes_duringTTL=[spikes_duringTTL,round(all_spk_ts(inds_duringTTL) * fs)']; %spike times of interest, in samples.
    clusterIDs_duringTTL=[clusterIDs_duringTTL, all_spk_IDs(inds_duringTTL)'];
end

wfparams.spikeTimes=spikes_duringTTL';  %
wfparams.spikeClusters = clusterIDs_duringTTL';  %IDs of each spike, when all spikes are listed in one vector

wf_laser=getWaveForms_with_std_and_bad(wfparams);   %extract the waveforms

%plot, compare against the baseline waveform extracted for the spikestruct.
% NB this will need to be redone for units that only really fired during
% the laser.
nplots=ceil(sqrt(length(wf_laser.unitIDs)));
wf_fig=figure('Color', 'w');

for g=1:length(wf_laser.unitIDs)
  
    this_clust=wf_laser.unitIDs(g);
    clust_ID=find(spikeStruct.cids==this_clust);
    c_chan=spikeStruct.c_channel(clust_ID);
    
    wf_mean_laser=wf_laser.waveFormsMean(g, c_chan, :);
    n_spks_used=sum(~isnan(wf_laser.spikeTimeKeeps(g,:)));
    wf_sem_laser=wf_laser.waveFormsSTD(g, c_chan, :)./sqrt(n_spks_used);
    
    wf_mean_all=spikeStruct.av_waveform{clust_ID};
    wf_sem_all=spikeStruct.std_waveform{clust_ID}/sqrt(spikeStruct.nWFs_extracted(clust_ID));
    
    wave_time=1000/fs *(1:length(wf_mean_all)); %convert to milliseconds
    
    subplot(nplots, nplots, g);
    shadedErrorBar(wave_time, wf_mean_all, wf_sem_all, 'b', 1);
    hold on
    shadedErrorBar(wave_time, wf_mean_laser, wf_sem_laser, 'r', 1);
    xlabel('Time (ms)');
    ylabel('\muV');
    title(['Cluster #' num2str(clust_ID)], 'FontWeight', 'normal');
end
   
end

