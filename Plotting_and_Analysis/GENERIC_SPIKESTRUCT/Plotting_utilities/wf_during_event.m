function [wf_fig] = wf_during_event(spikeStruct,laser_events, duration, fn)
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

wfparamsNOT=wfparams; %for extracting waveforms NOT during events

%make a table of times to consider,e.g during TTLs
spikes_duringTTL=[];
clusterIDs_duringTTL=[];  
all_spk_ts=spikeStruct.st;
all_spk_IDs=spikeStruct.clu;
duration_s=duration/1000;
all_spk_tsNOT=all_spk_ts; %we will cut out any spikes during laser from these
all_spk_IDsNOT=all_spk_IDs;
all_inds_duringTTLs=[];

for p=1:length(laser_events)   %get a list of spikes / cluster IDs which happened during laser.
    event_time=laser_events(p);
    inds_duringTTL=find(all_spk_ts>event_time & all_spk_ts< (event_time+duration_s) ); 
    spikes_duringTTL=[spikes_duringTTL,round(all_spk_ts(inds_duringTTL) * fs)']; %spike times of interest, in samples.
    clusterIDs_duringTTL=[clusterIDs_duringTTL, all_spk_IDs(inds_duringTTL)'];
    
    %spike times in the window of interest...
    all_inds_duringTTLs=[all_inds_duringTTLs; inds_duringTTL];     
end

fprintf('\n Extracting waveforms during specified events...')
wfparams.spikeTimes=spikes_duringTTL';  %
wfparams.spikeClusters = clusterIDs_duringTTL';  %IDs of each spike, when all spikes are listed in one vector
wf_event=getWaveForms_with_std_and_bad(wfparams);   %extract the waveforms

fprintf('\n Extracting waveforms at all other times...')
all_spk_tsNOT(all_inds_duringTTLs)=[];
all_spk_IDsNOT(all_inds_duringTTLs)=[];
wfparamsNOT.spikeTimes=round(all_spk_tsNOT * fs)
wfparamsNOT.spikeClusters=all_spk_IDsNOT;
wf_NOTevent=getWaveForms_with_std_and_bad(wfparamsNOT);   %extract the waveforms

%plot, compare against the baseline waveform extracted for the spikestruct.
% NB this will need to be redone for units that only really fired during
% the laser.
nplots=ceil(sqrt(length(wf_event.unitIDs)));
wf_fig=figure('Color', 'w');

for g=1:length(wf_event.unitIDs)
  
    this_clust=wf_event.unitIDs(g);
    clust_ID=find(spikeStruct.cids==this_clust);
    c_chan=spikeStruct.c_channel(clust_ID);
    
    wf_mean_event=wf_event.waveFormsMean(g, c_chan, :);
    n_spks_used=sum(~isnan(wf_event.spikeTimeKeeps(g,:)));
    wf_sem_event=wf_event.waveFormsSTD(g, c_chan, :)./sqrt(n_spks_used);
    
    this_clust_NOT=find(wf_NOTevent.unitIDs==this_clust)  
    wf_mean_NOTevent=wf_NOTevent.waveFormsMean(this_clust_NOT, c_chan, :);
    n_spks_usedNOT=sum(~isnan(wf_NOTevent.spikeTimeKeeps(this_clust_NOT,:)));
    wf_sem_NOTevent=wf_NOTevent.waveFormsSTD(this_clust_NOT, c_chan, :)./sqrt(n_spks_usedNOT);
       
    wave_time=1000/fs *(1:length(wf_mean_event)); %convert to milliseconds
    
    subplot(nplots, nplots, g);
    shadedErrorBar(wave_time, wf_mean_NOTevent, wf_sem_NOTevent, 'b', 1);
    hold on
    shadedErrorBar(wave_time, wf_mean_event, wf_sem_event, 'r', 1);
    xlabel('Time (ms)');
    ylabel('\muV');
    xlim([0, 2.7])
    title(['Cluster #' num2str(clust_ID)], 'FontWeight', 'normal');
     
    if g==length(wf_event.unitIDs)
        text(4,0, {'Red: waveform during laser' 'Blue: waveforms at other times'})
    end
end
   
end

