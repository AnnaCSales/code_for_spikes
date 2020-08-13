function wf = getWaveForms_withstd_and_bad(gwfparams)


% function wf = getWaveForms(gwfparams)
% TAKES ACCOUNT OF BAD (DISCONNECTED CHANNELS) - call the function with the
% additional field nBad = the number of bad channels not used in clustering
% Extracts individual spike waveforms from the raw datafile, for multiple
% clusters. Returns the waveforms and their means within clusters.
%
% This version also returns a standard deviation on the mean waveform.
% Contributed by C. Schoonover and A. Fink
%
% % EXAMPLE INPUT

% gwfparams.dataDir = '/path/to/data/';    % KiloSort/Phy output folder
% gwfparams.fileName = 'data.dat';         % .dat file containing the raw 
% gwfparams.dataType = 'int16';            % Data type of .dat file (this should be BP filtered)
% gwfparams.nCh = 32;                      % Number of channels that were streamed to disk in .dat file
% gwfparams.wfWin = [-40 41];              % Number of samples before and after spiketime to include in waveform
% gwfparams.nWf = 2000;                    % Number of waveforms per unit to pull out
% gwfparams.spikeTimes =    [2,3,5,7,8,9]; % Vector of cluster spike times ***(in samples)*** same length as .spikeClusters
% gwfparams.spikeClusters = [1,2,1,1,1,2]; % Vector of cluster IDs (Phy nomenclature)   same length as .spikeTimes
% ADDED:
% gwfparams.nBad  This is the number of disconnected channels.


% % OUTPUT
% wf.unitIDs                               % [nClu,1]            List of cluster IDs; defines order used in all wf.* variables
% wf.spikeTimeKeeps                        % [nClu,nWf]          Which spike times were used for the waveforms
% wf.waveForms                             % [nClu,nWf,nCh,nSWf] Individual waveforms
% wf.waveFormsMean                         % [nClu,nCh,nSWf]     Average of all waveforms (per channel)
%                                          % nClu: number of different clusters in .spikeClusters
%                                          % nSWf: number of samples per waveform
%
% % USAGE
% wf = getWaveForms(gwfparams);

%outfrom KS will be 32-nBad channels
% Load .dat and KiloSort/Phy output
fileName = fullfile(gwfparams.dataDir,gwfparams.fileName);     
filenamestruct = dir(fileName);
dataTypeNBytes = numel(typecast(cast(0, gwfparams.dataType), 'uint8')); % determine number of bytes per sample
nSamp = filenamestruct.bytes/(gwfparams.nCh*dataTypeNBytes);  % Number of samples per channel
wfNSamples = length(gwfparams.wfWin(1):gwfparams.wfWin(end));
mmf = memmapfile(fileName, 'Format', {gwfparams.dataType [gwfparams.nCh nSamp] 'x'});
chMap = readNPY([gwfparams.dataDir 'channel_map.npy'])+1;               % Order in which data was streamed to disk; must be 1-indexed for Matlab

% Read spike time-centered waveforms
unitIDs = unique(gwfparams.spikeClusters);
numUnits = size(unitIDs,1);
spikeTimeKeeps = nan(numUnits,gwfparams.nWf);
waveForms = nan(numUnits,gwfparams.nWf,gwfparams.nCh-gwfparams.nBad,wfNSamples);  %Amended to take account of disconnected channels.
waveFormsMean = nan(numUnits,gwfparams.nCh-gwfparams.nBad,wfNSamples);
waveFormsSTD = nan(numUnits,gwfparams.nCh-gwfparams.nBad,wfNSamples);
for curUnitInd=1:numUnits
    curUnitID = unitIDs(curUnitInd);  %specify which cluster we're working with
    curSpikeTimes = gwfparams.spikeTimes(gwfparams.spikeClusters==curUnitID); %extract all spike times for current cluster
    curUnitnSpikes = size(curSpikeTimes,1);  %number of spikes in this cluster
    spikeTimesRP = curSpikeTimes(randperm(curUnitnSpikes));   %randomly mix up the spike times (to allow us to select 2000 of them at random)
    spikeTimeKeeps(curUnitInd,1:min([gwfparams.nWf curUnitnSpikes])) = sort(spikeTimesRP(1:min([gwfparams.nWf curUnitnSpikes])));  %the minimum allows us to take everything there is, if the total num_spikes <2000

    for curSpikeTime = 1:min([gwfparams.nWf curUnitnSpikes])
        tmpWf = mmf.Data.x(1:gwfparams.nCh,spikeTimeKeeps(curUnitInd,curSpikeTime)+gwfparams.wfWin(1):spikeTimeKeeps(curUnitInd,curSpikeTime)+gwfparams.wfWin(end));
        waveForms(curUnitInd,curSpikeTime,:,:) = tmpWf(chMap,:);
    end
    waveFormsMean(curUnitInd,:,:) = squeeze(nanmean(waveForms(curUnitInd,:,:,:)));
    waveFormsSTD(curUnitInd,:,:) = squeeze(nanstd(waveForms(curUnitInd,:,:,:)));
    disp(['Completed ' int2str(curUnitInd) ' units of ' int2str(numUnits) '.']);
end

% Package in wf struct
wf.unitIDs = unitIDs;
wf.spikeTimeKeeps = spikeTimeKeeps;
wf.waveForms = waveForms;
wf.waveFormsMean = waveFormsMean;
wf.waveFormsSTD = waveFormsSTD;

end