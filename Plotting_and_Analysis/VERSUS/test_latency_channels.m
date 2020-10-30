%Quick script to test the idea of 'ranking' channels, based on whether they
%are displaying constant latency responses. Idea is that you could then
%pull out this channel for use with the latency tracker, 
%Uses pre-recorded data.

%path to data
% datapath='C:\Versus data\290920\Rec5\firstTwenty\';
datapath='C:\Versus data\281020\2020-10-28_16-13-14\CAR\';

%path to channel map
chanmap_fn='C:\MATLABcode\code_for_spikes\DesktopClustering\chanMapPoly3_25s.mat';
chanMap_struct=load(chanmap_fn);
chanMap=chanMap_struct.chanMap;
chanMap=fliplr(chanMap); %want it to start at the top and go down
%retrieve the TTLs
% TTLs = returnTTLs_not_from_zero(datapath);
fsts=footshock_TTLs(footshock_TTLs>1640);
%process footshock TTLs - these are in pairs (start, end) but we only want
%the one marking the onset.
% fsts=TTLs.digital{7};
% fsts(2:2:end)=[];

% Pull out sets delivered at difference frequencues

inds2=find( diff(fsts)>1.9 & diff(fsts) <2.1)+1;
inds2=[inds2(1)-1; inds2];
fsts_2Hz=fsts(inds2);


inds025=find( diff(fsts)>3.9 & diff(fsts) <4.1)+1;
inds025=[inds025(1)-1; inds025];
fsts_025Hz=fsts(inds025);

% now read the first 10s of a demo file to find the 
cnt_path='C:\Versus data\281020\2020-10-28_16-13-14\CAR\'; %path to the full list of continuous files.
 
% Pull out the postfix and prefix of then names of the continuous files
file_list=dir(cnt_path);
chk4cont = regexp({file_list.name}, '.continuous', 'once');
cont_matches=find(~cellfun(@isempty, chk4cont));

if length(cont_matches)>0
    
    fprintf('\n Continuous data...reading example file....\n')
    ex_file=file_list(cont_matches(1)).name;
    %only read the first few samples - this will give us the start time,
    %and metadata
    [ex_data, ts_start, info] = load_open_ephys_data([datapath, ex_file],'Indices', 1);   
    name_parts=strsplit(ex_file, '_');
    file_prefix=[name_parts{1} '_'];%the bits around 'CH1', i.e. 113_CH1_2, pre is 113_, post is _2
    
    if numel(name_parts)>2
        post_=strsplit(name_parts{3}, '.');
      
        if numel(post_)>1
           file_postfix=['_' post_{1}];
        else
           file_postfix=[];
        end
    else
          file_postfix=[];
    end    
end

%% Now go through the TTLs and store a window of data for each channel.
fs=30000;
window_length=400; %duration of window in ms
window_sample_length=(fs*0.001*window_length);  %# samples in window
nchans=32;  %' channels
tbase_window=1/fs*(1:window_sample_length);  %a time vector for the window

%specify which TTLs to use:
ttls_of_interest=fsts_025Hz;
nTTLs=10;  %how many to use for averaging
ttl_list=randperm(length(ttls_of_interest))  %pull out the required number of TTLs at random
ttl_list(nTTLs+1:end)=[];

ttls_to_use=ttls_of_interest(ttl_list); %randomly picked list.

%this variable will store the extracted data:
windowed_data=zeros(nchans,window_sample_length,nTTLs);

%now extract the data
for ttl=1:nTTLs
   %keep track of where we are
%    if ~mod(ttl, 2)
       fprintf('Working on TTL %d \n', ttl);
%    end
   
   this_TTL=ttls_to_use(ttl);
   TTL_in_samples=(this_TTL-ts_start) * fs; %time since start of recording * sampling rate.
   % work out the indices to extract from each file: take window from the
   % time of the TTL to 400ms later

   inds_to_include=fix([TTL_in_samples : (TTL_in_samples + window_sample_length -1) ]);
   
   for chan=1:32
      chan_name=[cnt_path file_prefix 'CH' num2str(chan) file_postfix '.continuous'];
      [data_win, ts_win,  ~] = load_open_ephys_data(chan_name,'Indices', inds_to_include); %read in first second only
      windowed_data(chan,:, ttl)=data_win;
      spike_times=spikes_by_threshold(data_win, fs, nstds);      
   end
end
%% Now extract spike crossings

skip_time=10;  %time window in ms to skip at the start of each event - for excluding artefacts.
binwin=0.002;  %for binning spikes
skip_bins=skip_time*0.001/(binwin);

binedges=(skip_time*0.001):binwin:0.3; %bin edges used for each window.
bin_cents=binedges;
bin_cents(end)=[];

%set up some labels for plotting
start_ind=find(bin_cents==0.05)
gap_ind=0.05/binwin
labels_x=[bin_cents(1), bin_cents(start_ind:gap_ind:end)];
xtick_ind=[1, [start_ind:gap_ind:length(bin_cents)]];

hist_data=zeros(nchans, numel(bin_cents), nTTLs);
nstds=3; % threshold crossing - for spike extraction

for ttl=1:nTTLs
   
   for chan=1:32
     
      data_=windowed_data(chan,:, ttl);
      spike_times=spikes_by_threshold(data_, fs, nstds);
      [bindata, ~]=histcounts(spike_times, binedges);
      hist_data(chan,:,ttl)=bindata;
      
   end
end


meandata=mean(hist_data,3);
overall_mean=mean(reshape(meandata, numel(meandata), []));
z_meandata=zscore(meandata, 0,2)

% figure('Color', 'w')
% % bar(max(z_meandata(:,8:end), [],2))
% xticks(1:32)
% set(gca,'fontsize',8)
% xlabel('OEP channel number')
% ylabel('Max value of zscored mean binned #spikes')

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
%% Try again binning voltage, no spike extraction
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

volt_bin_data=zeros(nchans,(window_sample_length/n_samples_per_bin)-skip_bins+1, nTTLs);
%  volt_bin_data=[];

for ttl=1:nTTLs
   
   for chan=1:32
     
      data_=abs(windowed_data(chan,:, ttl));
      bin_data=reshape(data_, n_samples_per_bin, [])';
      mean_bin_data=mean(bin_data(skip_bins:end, :),2);    
      volt_bin_data(chan,:,ttl)=mean_bin_data;
      
   end
end

meandata=mean(volt_bin_data,3);
z_meandata=zscore(meandata, 0,2)

% figure('Color', 'w')
% for r=1:32
% this_chan=chanMap(r);
% hold on
% plot(( z_meandata(r,:))+20*r)
% end
% 
%  yticks(14:20:20*nchans)
% yticklabels(fliplr(chanMap))
% ylabel('OEP channel #')
% % ylim([0,9.1*nchans])
% title(['Mean binned voltage over ' num2str(nTTLs) ' trials'], 'FontWeight', 'normal')
% 

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
caxis([0,7])
title(['Average reponse to ' num2str(nTTLs) ' trials'])