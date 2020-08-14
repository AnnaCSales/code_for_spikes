function [sdf, tbase] = unitSDF_all(spikeStruct, cellIDs, time, gaussSDtime)
 %Merges spikes from listed clusters into one big spike train
% Returns a smooth spike density function by convolving with a gaussian
% kernel. Expects: spikeStruct, the cluster number of cell to plot, the
% time period ( a two element vector of the form [min max]), and the width
% of the gaussian kernel in milliseconds.
all_ts=[];
for cl=1:length(cellIDs) %merge  
    ts_= spikeStruct.timesSorted{cellIDs(cl)};
    all_ts=[all_ts; ts_];
end

all_ts=sort(all_ts, 'ascend');
binsize= 0.001; 
tbin_edges= time(1):binsize:time(2);  %look at FR properties in baseline only

spk_count_unit = histc(all_ts,tbin_edges);
spk_count_unit = spk_count_unit(1:end-1);
gauss_window = 1./binsize; % 1 second window
SDtime=gaussSDtime/1000; %convert width to seconds.
gauss_SD = SDtime./binsize; % 50ms SD on the gaussian
gk = gausskernel(gauss_window,gauss_SD); 
gk = gk./binsize; % normalize by binsize

sdf = conv2(spk_count_unit,gk,'same'); % convolve with gaussian window
tbase = tbin_edges(1:end-1)+binsize/2;

sdf(end)=[]; %end up with one too many due to the way the edges of bins are constructed.
tbase(end)=[];
end

