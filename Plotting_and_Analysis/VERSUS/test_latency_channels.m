% This is a demo script to test the idea of using raw data to identify channels 
% showing spiking responses at constant latency. The idea is to explore 
% processing options that could be implemented in real time on OEP. This would 
% help to identify channels of interest during an experiment (i.e. those
% that are most likely to show C fibre activity. This version uses
% pre-recorded data to test different approaches

% Anna Sales, November 2020

%% Specify path to prerecorded data
datapath='C:\Versus data\281020\2020-10-28_16-13-14\CAR\';

% specify a path to a channel map file (script expects this in kilosort
% format, i.e a struct with a field called 'chanMap', which lists OEP channels
% from bottom of probe upwards
chanmap_fn='C:\MATLABcode\code_for_spikes\DesktopClustering\chanMapPoly3_25s.mat';
chanMap_struct=load(chanmap_fn);
chanMap=chanMap_struct.chanMap;
nchans=length(chanMap);  %' number of channels

%% Retrieve the TTLs marking electrical stimulation (footshock timestamps - 'fsts')

% either extract from the raw data:

% TTLs = returnTTLs_not_from_zero(datapath);
%process footshock TTLs - these are in pairs (start, end) but we only want
%the one marking the onset:
% fsts=TTLs.digital{7};
% fsts(2:2:end)=[];


% Or use preprocessed file:
fsts=footshock_TTLs;

% Pull out sets of shocks delivered at difference frequencues

inds2=find( diff(fsts)>0.48 & diff(fsts) <0.51)+1;
inds2=[inds2(1)-1; inds2];
fsts_2Hz=fsts(inds2);  % 2Hz stims


inds025=find( diff(fsts)>3.9 & diff(fsts) <4.1)+1;
inds025=[inds025(1)-1; inds025];
fsts_025Hz=fsts(inds025);  %0.25Hz stims

% Start to read in data - work out the format of the continuous file names
% in the folder.
[file_prefix, file_postfix]=returnOEPformat(datapath);

%% Now go through the TTLs and store a window of data for each channel.
fs=30000;  %sampling rate
window_length=300; %duration of window after each shock to keep, in ms
window_sample_length=(fs*0.001*window_length);  % number of samples in the window

tbase_window=1/fs*(1:window_sample_length);  %a time vector for the window

%specify which TTLs to use in the analysis:
ttls_of_interest=fsts_2Hz;
nTTLs=150;  %how many of them to use 
% ttl_list=randperm(length(ttls_of_interest))  %option - pull out the required number of TTLs at random from the set.
ttl_list=1:nTTLs;
if nTTLs<length(ttls_of_interest)  %
    ttl_list(nTTLs+1:end)=[];
    ttls_to_use=ttls_of_interest(ttl_list); %randomly picked list.
else
    fprintf('\n Number of TTLs higher than number in set of interest - using all TTLs \n')
    ttls_to_use=ttls_of_interest;
end

%this variable will store the extracted data:
windowed_data=zeros(nchans,window_sample_length,nTTLs);

%now extract the data. NB this is a slow process with >~10 TTLs.
for chan=1:32

   fprintf('Extracting TTL windows on channel %d \n', chan);

   chan_name=[datapath file_prefix 'CH' num2str(chan) file_postfix '.continuous'];
   [data_full, ts_full,  ~] = load_open_ephys_data(chan_name); %read in entire chan
 
   for ttl=1:nTTLs   
      this_TTL=ttls_to_use(ttl);
      ttl_win_times=[this_TTL, this_TTL+(0.001*window_length)];
      inds_in_window=find(ts_full>=ttl_win_times(1) & ts_full<=ttl_win_times(2));
      data_win = data_full(inds_in_window);
      data_win(window_sample_length+1:end)=[]; %sometimes ends up with one sample too many
      try
          windowed_data(chan,:, ttl)=data_win;
      catch
          fprintf('\n Error windowing data \n')
      end
   end
end
%% (1) Averaging over trials, look for 'hotspots' on the probe where constant latency responses are occuring.
% First try to do this by extracting spikes by simple threshold crossing,
% and counting spikes in bins.

skip_time=10;  %time window in ms to skip at the start of each event - for excluding artefacts.
binwin=2;  %bn width in m seconds, for binning spikes
binwin=0.001*binwin;
skip_bins=skip_time*0.001/(binwin);

binedges=(skip_time*0.001):binwin:0.3; %bin edges used for each window.
bin_cents=binedges;
bin_cents(end)=[];

%set up some labels for plotting
start_ind=find(bin_cents==0.05);
gap_ind=0.05/binwin;
labels_x=[bin_cents(1), bin_cents(start_ind:gap_ind:end)];
xtick_ind=[1, [start_ind:gap_ind:length(bin_cents)]];

%hist_data is a variable which will store binned spike counts for each TTL on
%each channel:
hist_data=zeros(nchans, numel(bin_cents), nTTLs);

nstds=3; % sets the threshold for spike extraction in terms of number of stds from mean.

for ttl=1:nTTLs   
   for chan=1:32
     
      data_=windowed_data(chan,:, ttl);
      spike_times=spikes_by_threshold(data_, fs, nstds);
      [bindata, ~]=histcounts(spike_times, binedges);
      hist_data(chan,:,ttl)=bindata;
      
   end
end

% Average over trials. For each channel, the average should begin to
% reflect level of background noise over multiple trials.
meandata=mean(hist_data,3);

% Now z-score, to find above average bins for a given channel, and to make 
% comparisons between channels easy to do.
z_meandata=zscore(meandata, 0,2)


% Plot the data across channels, to see if we can pull out channels of
% interest.
figure('Color', 'w')
for r=1:32 %start plotting at the lowest channel on the probe, work upwards
    this_chan=chanMap(r);
    hold on
    plot(bin_cents,( meandata(this_chan,:))+r)
end

yticks(1:32)
yticklabels(chanMap)
ylabel('OEP channel #')
ylim([0,33])
title(['Mean binned spike count over ' num2str(nTTLs) ' trials'], 'FontWeight', 'normal')

figure('Color', 'w')
imagesc(meandata(fliplr(chanMap),:))
xticks(xtick_ind)
xticklabels(string(labels_x))
xlabel('Time (ms)')
ylabel('Channel')
yticks(1:32)
yticklabels(fliplr(chanMap))
title(['Average reponse to ' num2str(nTTLs) ' trials'])
%% Try the same approach, but this time just use the average raw voltage in each time bin
%  and don't bother binning and counting spikes. Over multiple trials,
%  peaks in voltage at given times should become apparent.

binwin2=0.001;  %binning window for voltage.
skip_time=10;  %time window in ms to skip at the start of each event - for excluding artefacts.
skip_bins=round((skip_time*0.001)/binwin2);  %number of bins to skip.
binedges2=(skip_time*0.001):binwin2:(0.001*window_length); %bin edges used for each window.
bin_cents2=binedges2;
bin_cents2(end)=[];
n_samples_per_bin=binwin2*fs;


%set up some labels for plotting
start_ind2=find(bin_cents2==0.05)
gap_ind2=0.05/binwin2
labels_x=[bin_cents2(1), bin_cents2(start_ind2:gap_ind2:end)];
xtick_ind=[1, [start_ind2:gap_ind2:length(bin_cents2)]];

% 'volt_bin_data' will hold the binned voltage data.
volt_bin_data=zeros(nchans,(window_sample_length/n_samples_per_bin)-skip_bins+1, nTTLs);

for ttl=1:nTTLs 
   for chan=1:32
     
      data_=abs(windowed_data(chan,:, ttl));
      bin_data=reshape(data_, n_samples_per_bin, [])';
      mean_bin_data=mean(bin_data(skip_bins:end, :),2);    
      volt_bin_data(chan,:,ttl)=mean_bin_data;
      
   end
end

% As before, average over trials, then z-score
meandata=mean(volt_bin_data,3);
z_meandata=zscore(meandata, 0,2)

figure('Color', 'w')
for r=1:32
this_chan=chanMap(r);
hold on
plot(( z_meandata(r,:))+20*r)
end

 yticks(14:20:20*nchans)
yticklabels(fliplr(chanMap))
ylabel('OEP channel #')
% ylim([0,9.1*nchans])
title(['Mean binned voltage over ' num2str(nTTLs) ' trials'], 'FontWeight', 'normal')


figure('Color', 'w')
bar(max(z_meandata(:,8:end), [],2))
xticks(1:32)
set(gca,'fontsize',8)
xlabel('OEP channel number')
ylabel('Max value of zscored mean binned voltage')


figure('Color', 'w')
imagesc(z_meandata(fliplr(chanMap),:))
xticks(xtick_ind)
 xticklabels(string(labels_x))
xlabel('Time (ms)')
aa=gca
xlim([0.5, 0.3/binwin2])
ylabel('Channel')
yticks(1:32)
yticklabels(fliplr(chanMap))
caxis([0,6])
title(['Average reponse to ' num2str(nTTLs) ' trials'])

%% Now try producing heatmaps, but flattening over channels. 
%  This will produce plots of latency vs TTL number - can look for latency
%  changes happening anywhere on the probe.

%  Try first with spike extraction, as above.
nEvents=50; %for plotting a smaller number of events than we extracted

if ~exist('nTTLs')
    nTTLs=size(windowed_data,3)
end

%extract spike crossings
skip_time=10;  %time window in ms to skip at the start of each event - for excluding artefacts.
binwin=0.001;  %for binning spikes
skip_bins=skip_time*0.001/(binwin);
binedges=(skip_time*0.001):binwin:0.3; %bin edges used for each window.
bin_cents=binedges;
bin_cents(end)=[];

%set up some labels for plotting
start_ind=find(bin_cents==0.05)
gap_ind=0.05/binwin
labels_x=[bin_cents(1), bin_cents(start_ind:gap_ind:end)];
xtick_ind=[1, [start_ind:gap_ind:length(bin_cents)]];

hist_data=zeros(nEvents, numel(bin_cents),nchans ); %for storing binned spike counts
nstds=3.4; % threshold crossing - for spike extraction

amp_data=[];  %additional variable to hold amplitude info, as well as whether a spike has occured.

for ttl=1:nEvents
   
   for chan=1:32
      binamps=zeros(size(bin_cents));
      data_=windowed_data(chan,:, ttl);
      [spike_times, spike_amps]=spikes_by_threshold(data_, fs, nstds);
      [bindata, ~, binID]=histcounts(spike_times, binedges);
      hist_data(ttl,:,chan)=bindata;
      if sum(binID)
          spike_amps(binID==0)=[]; %don't keep the ones that happened before skipcount
          binID(binID==0)=[];
          binamps(binID)=(spike_amps);  %amplitude associated with each bin, if there is a spike there.
      end     
      amp_data(ttl,:,chan)=binamps;  %instead of a '1' when there's a spike, this has the amplitude.
   end
end

%Plot heatmaps for each quadrant of the probes. This will enable
%determination of roughly where the activity is

spikes_quad=[];
quad1=chanMap(1:8);
spikes_quad(:,:,1)=sum( hist_data(:,:,quad1), 3);

quad2=chanMap(9:16);
spikes_quad(:,:,2)=sum( hist_data(:,:,quad2), 3);

quad3=chanMap(17:24);
spikes_quad(:,:,3)=sum( hist_data(:,:,quad3), 3);

quad4=chanMap(25:32);
spikes_quad(:,:,4)=sum( hist_data(:,:,quad4), 3);

title_labels={'Quad 1 (bottom)', 'Quad 2', 'Quad 3', 'Quad 4 (top)'};
figure('Color', 'w')
for q=1:4
    subplot(2,2,q)
    imagesc(flipud(spikes_quad(:,:, q)));
    xticks(xtick_ind)
    xticklabels(string(labels_x))
    xlabel('Time (ms)')
    yticks(1:5:nEvents)
    yticklabels(string(nEvents:-5:1))
    ylabel('Trial')
    title(title_labels{q})
    xlim([1,100])
end

%% Finally do the same thing, but with the raw voltage (no spike extraction)

if ~exist('nTTLs')
    nTTLs=size(windowed_data,3)
end

binwin2=0.001;
skip_time=10;  %time window in ms to skip at the start of each event - for excluding artefacts.
skip_bins=round((skip_time*0.001)/binwin2);  %number of bins to skip.

binedges2=(skip_time*0.001):binwin2:(0.001*window_length); %bin edges used for each window.
bin_cents2=binedges2;
bin_cents2(end)=[];


%set up some labels for plotting
start_ind2=find(bin_cents2==0.05)
gap_ind2=0.05/binwin2
labels_x=[bin_cents2(1), bin_cents2(start_ind2:gap_ind2:end)];
xtick_ind=[1, [start_ind2:gap_ind2:length(bin_cents2)]];

n_samples_per_bin=binwin2*fs;

volt_bin_data=zeros( nTTLs,(window_sample_length/n_samples_per_bin)-skip_bins+1,nchans);
%  volt_bin_data=[];

for ttl=1:nTTLs
   
   for chan=1:32
     
      data_=abs(windowed_data(chan,:, ttl));
      bin_data=reshape(data_, n_samples_per_bin, [])';
      mean_bin_data=mean(bin_data(skip_bins:end, :),2);    
      volt_bin_data(ttl,:,chan)=mean_bin_data;
      
   end
end

spikes_quad=[];
quad1=chanMap(1:8);
spikes_quad(:,:,1)=sum( volt_bin_data(:,:,quad1), 3);

quad2=chanMap(9:16);
spikes_quad(:,:,2)=sum( volt_bin_data(:,:,quad2), 3);

quad3=chanMap(17:24);
spikes_quad(:,:,3)=sum( volt_bin_data(:,:,quad3), 3);

quad4=chanMap(25:32);
spikes_quad(:,:,4)=sum( volt_bin_data(:,:,quad4), 3);

title_labels={'Quad 1 (bottom)', 'Quad 2', 'Quad 3', 'Quad 4 (top)'};
figure('Color', 'w')
for q=1:4
    subplot(2,2,q)
    imagesc(flipud(spikes_quad(:,:, q)));
    xticks(xtick_ind)
    xticklabels(string(labels_x))
    xlabel('Time (ms)')
    yticks(1:1:nTTLs)
    yticklabels(string(nTTLs:-1:1))
    ylabel('Trial')
    title(title_labels{q})
     xlim([0,200])
    caxis([30,80])
end