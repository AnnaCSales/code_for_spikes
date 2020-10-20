%Quick script to test the idea of 'ranking' channels, based on whether they
%are displaying constant latency responses. Idea is that you could then
%pull out this channel for use with the latency tracker, 
%Uses pre-recorded data.

%path to data
datapath='C:\Versus data\290920\Rec5\firstTwenty\';

%retrieve the TTLs
TTLs = returnTTLs_not_from_zero(datapath);

%process footshock TTLs - these are in pairs (start, end) but we only want
%the one marking the onset.
fsts=TTLs.digital{3};
fsts(2:2:end)=[];

% Pull out sets delivered at difference frequencues
fsts_2Hz=fsts;
exc2Hz=find(diff(fsts_2Hz)<0.49 | diff(fsts_2Hz)>0.51); 
fsts_2Hz(exc2Hz+1)=[];

fsts_025Hz=fsts;
exc025=find(diff(fsts_025Hz)<4| diff(fsts_025Hz) >4.5 );
fsts_025Hz(exc025+1)=[];

% now read the first 10s of a demo file to find the 
cnt_path='C:\Versus data\290920\Rec5\'; %path to the full list of continuous files.

% Pull out the postfix and prefix of then names of the continuous files
file_list=dir(cnt_path);
chk4cont = regexp({file_list.name}, '.continuous', 'once');
cont_matches=find(~cellfun(@isempty, chk4cont));

if length(cont_matches)>0
    
    fprintf('\n Continuous data...reading example file....\n')
    ex_file=file_list(cont_matches(1)).name;
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

window_length=400; %duration of window in ms
window_sample_length=(fs*0.001*window_length);  %# samples in window
nchans=32;  %' channels
tbase_window=1/fs*(1:window_sample_length);  %a time vector for the window

%this variable will store the extracted data:
windowed_data=zeros(nchans,window_sample_length,nTTLs);

nstds=2.5; % threshold crossing - for spike extraction

%specify which TTLs to use:
ttls_to_use=fsts_025Hz;
nTTLs=20;  %how many to use for averaging
ttl_list=randperm(length(ttls_to_use))  %pull out the required number of TTLs at random
ttl_list(nTTLs+1:end)=[];

%now extract the data
for ttl=1:nTTLs
   %keep track of where we are
   if ~mod(ttl, 2)
       fprintf('Working on TTL %d ', ttl);
   end
   
   this_TTL=ttls_to_use(ttl);
  
   % work out the indices to extract from each file: take window from the
   % time of the TTL to 400ms later

   inds_to_include=fix([this_TTL*fs : ((this_TTL*fs) + window_sample_length -1) ]);
   
   for chan=1:32
      chan_name=[cnt_path file_prefix 'CH' num2str(chan) file_postfix '.continuous'];
      [data_win, ts_win,  ~] = load_open_ephys_data(chan_name,'Indices', inds_to_include); %read in first second only
      windowed_data(chan,:, ttl)=data_win;
      spike_times=spikes_by_threshold(data_win, fs, nstds);      
   end
end
%% Now extract spike crossings

binwin=0.0005;  %for binning spikes
binedges=0.003:binwin:0.3; %there is an artefact lasting ~3ms
bin_cents=binedges;
bin_cents(end)=[];
hist_data=zeros(nchans, numel(bin_cents), nTTLs);

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

figure('Color', 'w')
% bar(max(z_meandata(:,8:end), [],2))
xticks(1:32)
set(gca,'fontsize',8)
xlabel('OEP channel number')
ylabel('Max value of zscored mean binned #spikes')

figure('Color', 'w')
for r=1:32
hold on
plot(bin_cents,( meandata(r,:))+r)
end

yticks(1:32)
yticklabels(1:32)
ylabel('OEP channel #')
ylim([0,33])
title(['Mean binned spike count over ' num2str(nTTLs) ' trials'], 'FontWeight', 'normal')
%% Try again binning voltage, no spike extraction
binwin2=0.001;

n_samples_per_bin=binwin2*fs;

volt_bin_data=zeros(nchans,window_sample_length/n_samples_per_bin, nTTLs);
for ttl=1:nTTLs
   
   for chan=1:32
     
      data_=abs(windowed_data(chan,:, ttl));
      bin_data=reshape(data_, n_samples_per_bin, [])';
      mean_bin_data=mean(bin_data,2);    
      volt_bin_data(chan,:,ttl)=mean_bin_data;
      
   end
end

meandata=mean(volt_bin_data,3);
z_meandata=zscore(meandata, 0,2)

figure('Color', 'w')
for r=1:32
hold on
plot(( meandata(r,:))+8*r)
end
 yticks(14:8:9*nchans)
 yticklabels(1:32)
ylabel('OEP channel #')
ylim([0,9.1*nchans])
title(['Mean binned voltage over ' num2str(nTTLs) ' trials'], 'FontWeight', 'normal')


figure('Color', 'w')
bar(max(z_meandata(:,8:end), [],2))
xticks(1:32)
set(gca,'fontsize',8)
xlabel('OEP channel number')
ylabel('Max value of zscored mean binned voltage')
