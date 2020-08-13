function [sdf, tbase] = unitSDF(spikeStruct, cellID, time, gaussSDtime)

% Returns a smooth spike density function by convolving with a gaussian
% kernel. Expects: spikeStruct, the cluster number of cell to plot, the
% time period ( a two element vector of the form [min max]), and the width
% of the gaussian kernel in milliseconds.
    ts_= spikeStruct.timesSorted{cellID};
    binsize= 0.001; 
    tbin_edges= time(1):binsize:time(2);  %look at FR properties in baseline only

    spk_count_unit = histc(ts_,tbin_edges);
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

