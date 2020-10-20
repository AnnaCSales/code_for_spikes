function [spike_times] = spikes_by_threshold(data, fs, nstds)
%Quick script to read in a window of data, cut out extreme values, extract spikes. Produces a
%vector of spike times. Follows Quiroga (2004, Neural Computation) for details of
%threshold and definitions of noise.

% 'nstds' is the number of stds from the mean
% 'fs' is the sampling rate, in Hz
% 'data' is the raw data 

tst=[1:length(data)]  * 1/fs; %timebase for the data
%%  Get rid of the worst extreme values
cutoff=1000;
cut_inds=find(data>cutoff | data<(-1*cutoff));  
data(cut_inds)=[];
tst(cut_inds)=[];

%% Now estimate the noise and set a threshold for spike detection

times_noise=nstds;
sd_est=median(abs(data)./0.6745) ;   %estimate of the noise, (following Quiroga(2004) )
threshold=-times_noise*sd_est ; %negative threshold
sd1=ones(numel(tst), 1).*threshold;  
sd2=-1.*sd1;

%%  Extract spikes.
%Look only at points above the threshold, find peaks 

inds_above=find(data<threshold);
data_above=data(inds_above);
[pks,locs] = findpeaks(-1.*data_above);  %find locations of minima

%now pull out the actual spikes.
peak_inds=inds_above(locs(:));  %index of spikes in full dataset
num_spikes=length(peak_inds);

%return the spike times.
spike_times=tst(inds_above(locs));  %spike timestamps.

end
