%Script takes output from phy and open ephys and produces a struct which contains
%clustered spike data, average waveforms, TTLs recorded in OEP.
%
%Updated July 2020 to automatically process and add TTL information from
%openephys in either binary or continuous format. If using binary data,
%make a folder in the directory containing the PHY output and put all TTL
%data in there (in binary format, OEP generates several files called
%'timestamps.npy' and if they are not in a separate folder, it may process
%the wrong one). Name this folder 'TTL'.

%Dependencies:
% returnTTLs - script to process digital TTLs.
% KSdir (OEP utility)  - which in turn requires loadParamsPy, readNPY
% (npy-matlab-master), from PHY or OEP GIT.
% Amended versions of getWaveForms / getWaveForms_withBad (PHY/OEP utilities)
% Add the folder containing dependencies to the PATH
% addpath(genpath('D:\Code\MATLAB\AnalysingPHYoutput'));
% addpath(genpath('D:\Code\MATLAB\Packages\KS_PHY_OEP_utilities\preprocessing\'));
% Provide the name of directory with the data in it. This dir should contain the binary file
% used for clustering, the PHY output, at least one OEP continuous file if using 
% continuous format / a folder called TTL if using binary. OEP metadata
% (events, messages files etc) should be in there too.
rootpath=[pwd '\'];
% Provide the name of the binary file used for clustering.
rawData='dataALL.bin';

% Provide the name of the folder where the TTL info is, if using binary
% (this field does not matter if you are using continuous format)
TTLfolder='TTL';

bl_start=0;  %define the baseline period, if there is one.
bl_end=1;

%% Extract TTLs. 
struct_file_pn='rootpath';
TTLs = returnTTLs(rootpath, struct_file_pn);

%% Pre-processing 
% run the OEP utility to create the basic struct
spikeStruct = loadKSdir(rootpath);  %ONLY CLUSTERS MARKED AS 'GOOD' in PHY ARE TAKEN INTO THE CODE

%add a few more useful fields:
spikeStruct.nclusts  = length(spikeStruct.cids);
spikeStruct.baseline_st=bl_start;
spikeStruct.baseline_end=bl_end;
spikeStruct.TTLs=TTLs; 
cgs=spikeStruct.cgs;
spikeStruct.baseline_st=bl_start;
spikeStruct.baseline_end=bl_end;
%% Extract more cluster info
fid = fopen('cluster_info.tsv');
cluster_info = textscan(fid, '%d %f %f %s %f %f %f %f %s %d %d', 'HeaderLines', 1);
fclose(fid);

fid2 = fopen('cluster_group.tsv');
C = textscan(fid2, '%s%s', 'HeaderLines', 1);
fclose(fid2);
chanFlag=0;
if ( length(cluster_info{1}) ~= length(C{1}))
    fprintf('\n *** Fix blanks in cluster_info before proceeding ****\n')
    fprintf('\nIf script is continued, will attempt to estimate centre channel\n')
    chanFlag=1;
else
    c_channels_all=cluster_info{6};
    contam_pc_all=cluster_info{3};
end

[~, groups] = readClusterGroupsCSV('cluster_group.tsv');
good_clusts=find(groups==2); %all the ones marked as 'good' (green) in PHY
c_channelsOEP=c_channels_all(good_clusts);  % this is given in terms of OEP chans
[~,c_channelPHY]=ismember(c_channelsOEP, spikeStruct.chanMap);

contam_pc=contam_pc_all(good_clusts);
spikeStruct.contam_pc=contam_pc;
%%
%Check for bad (disconnected) channels and flag if there are any - will affect waveform extraction later on
nChansBad=spikeStruct.n_channels_dat-length(spikeStruct.chanMap)   ;

% Get the sampling rate.
fs=spikeStruct.sample_rate;

%% Spike times
% extract times of cluster spikes and add them to the struct.
for iClu = 1:length(spikeStruct.cids)    
    spikeStruct.timesSorted{iClu} = spikeStruct.st(spikeStruct.clu==spikeStruct.cids(iClu)); % one cell per cluster
    [ro, co, ~]=find(spikeStruct.clu==spikeStruct.cids(iClu)) ; %find the index of all spikes from one cluster   
    spikeStruct.clu_templates{iClu}=spikeStruct.spikeTemplates(ro); %the template ID of each spike in the cluster (some will be a mixture)    
end
%% Waveform information

%Update the waveform for each cluster. Don't need to run this more than
%once 
nWF=2000;
if ~exist([rootpath , '\newWFs.mat'])  %if we haven't already done it, go get the new waveforms
    %parameters for the waveform getting widget.
    wfparams.dataDir=rootpath;
    wfparams.fileName=rawData;
    wfparams.nCh = 32;
    wfparams.nBad = nChansBad;
    wfparams.dataType='int16';
    wfparams.wfWin = [-40 41];
    wfparams.spikeClusters = spikeStruct.clu;
    wfparams.nWf = nWF; 
    wfparams.spikeTimes=round((spikeStruct.st) * fs);  %because we need this in SAMPLES not times. 
    
    %Sometimes the spike times are VERY slightly off (~10e-12)due to the errors that
    %creep in when blocks are split and merged etc- the 'round' makes sure that
    %everything in this vector is an integer.
    %lets check - run this if errors might be getting big
    [a,b]=find(mod(wfparams.spikeTimes(:), 1)~=0);
    if a
        for p=1:length(a)
            c(p)=round(wfparams.spikeTimes(a(p)))-wfparams.spikeTimes(a(p));   
        end
    else
        c=0;
    end
    fprintf(' \n Max error on spike timestamps: %f' , max(c)')
    fprintf(' \n Extracting waveforms.....');
    if nChansBad
        wf=getWaveForms_with_std_and_bad(wfparams);  %widget to extract a random sample of 2000 of the cluster wave forms and return all the individual traces and the average waveform on each channel
        fprintf('Disconnected channels detected, extracting waveforms on remaining channels only')
    else
        wf=getWaveForms_with_std_and_bad(wfparams); 
    end
    mean_wf=wf.waveFormsMean;
    std_wf=wf.waveFormsSTD;
    
    n_spks=[];
    for w=1:size(wf.spikeTimeKeeps,1) %go through unit, logging how many waveform examples were saved
      n_spks(w)=sum(~isnan(wf.spikeTimeKeeps(w,:)));
    end
    
    wf_info.meanWF=mean_wf;
    wf_info.stdWF=std_wf;
    wf_info.numSpksUsed=n_spks;    
    save([rootpath '\wf_info.mat'], 'wf_info');   %too big to save all of it, for now anyhow - just save the averages.
 

else
    wf_inf=load([rootpath, 'wf_info.mat']);
    mean_wf=wf_inf.wf_info.meanWF;
    std_wf= wf_inf.wf_info.stdWF;
    n_spks=wf_inf.wf_info.numSpksUsed;
end

spikeStruct.allchanWFs=mean_wf;  %this is the mean WF for all channels  
spikeStruct.allchanSTDs=std_wf;   %this is the STD on the WF for each channel.
spikeStruct.nWFs_extracted=n_spks; %this is the number of waveforms sampled to get the 

%% Channel information for each cluster
%extract the centre channel for the mean waveform on that centre channel for each cluster

for iUnit = 1:length(spikeStruct.cids)   
   this_wf = squeeze(mean_wf(iUnit, :, :)); %get rid of the singleton dim
   this_std= squeeze(std_wf(iUnit, :,:));
 
  if chanFlag %if PHY output has corrupted / not filled in correctly, try to estimate
   [~, c_channelest]=max(range(this_wf'));  %this will estimate the centre channel 
   spikeStruct.c_channel(iUnit)=c_channelest; %NB channels in Phy are also zero indexed so chan 1 in phy = chan 2 here.
   spikeStruct.av_waveform{iUnit}=this_wf(c_channelest, :);  %store the average waveform on the centre channel    
   spikeStruct.std_waveform{iUnit}=this_std(c_channelest, :);  %store the average waveform on the centre channel   
   fprintf('\n Estimating centre channel for unit \n' )
  else
    spikeStruct.c_channel(iUnit)=c_channelPHY(iUnit);
    spikeStruct.av_waveform{iUnit}=this_wf(c_channelPHY(iUnit), :);  %store the average waveform on the centre channel    
    spikeStruct.std_waveform{iUnit}=this_std(c_channelPHY(iUnit), :);  %store the average waveform on the centre channel   
  end
      
  
end

%% Useful vector of depths for plotting - clusters listed by depth order.
%work out a new order for displaying the spikes, based on their depth. If
%more than one group is centred around the same channel they will be
%at least be displayed next to each other.

[vals, indys]=sort(spikeStruct.c_channel);
plot_pos=[];
for iUnit = 1:1:length(spikeStruct.cids)
    
    [~,plot_pos(iUnit)]=find(indys==iUnit) ; % a plot position for each cluster, with lower (deeper) channels getting lower numbers
    
    %plot_pos is a vector n_clusts long with a position for each unit,
    %provided in the same order as usual.
end
spikeStruct.plot_pos=plot_pos;
%%
%pull out the limits of the recording - useful for all sorts of things.

for iUnit=1:1:spikeStruct.nclusts
    ts_= spikeStruct.timesSorted{iUnit};  %extract timestamps for this unit     
    %extract some more useful info about time range for plotting
    min_t2=min(ts_);
    max_t2=max(ts_);
    
    if iUnit==1 || min_t2<min_t
        min_t=min_t2;
    end
    
     if iUnit==1 || max_t2>max_t
        max_t=max_t2;
     end
end

spikeStruct.timeRange=[min_t, max_t]; %save the time range of the spikes.

%%
%Now save the spikeStruct
save( [rootpath '\spikeStruct.mat'], 'spikeStruct', '-v7.3');   
fprintf('\n Finished!')