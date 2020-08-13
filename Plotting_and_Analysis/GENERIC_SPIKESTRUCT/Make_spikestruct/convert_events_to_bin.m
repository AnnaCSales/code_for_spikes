% Script for converting an events file to a continuous binary file. Uses
% one continuous file recorded at the same time as a template, inserts 1s
% over a specified number of ms (set by 'TTL_length') 
% where there is a TTL, zero otherwise, in the resulting binary
% file

rootpath='D:\Versus\2807\2807oepformat\20200728\Rec1\CommonAvRef\'
[event_data, event_ts, event_info] = load_open_ephys_data([rootpath 'all_channels_2.events']);
[data, ts, info] = load_open_ephys_data([rootpath '115_CH1_2.continuous']);   
fs=info.header.sampleRate; % extract sampling rate

new_TTL_file=zeros(size(ts));

TTL_length=10; %length in ms of TTL pulse we want to write
num_samples_high=TTL_length * 0.001 * fs;

event_data_chan=2; %the event data we are interested in, i.e. event '1' or event '2' etc.
event_inds=find(event_data==event_data_chan); %index events we care about

for t=1:length(event_inds)
    this_event=event_ts(event_inds(t));    
    ts_ind=find(ts==this_event); %they should be exact matches as they were rec'd simultaneously 
    new_TTL_file(ts_ind:ts_ind+num_samples_high-1)=ones(num_samples_high,1);   
    
    if mod(t,50)==0 && t>1
        fprintf('\n Done %i of %i', t, numel(event_inds))
    end
end

 fileOut = fopen([rootpath, '/TTL_3.bin'],'w');
 fwrite(fileOut,new_TTL_file, 'int16');
 fclose(fileOut)
 
%sanity check
figure
plot(event_ts(event_inds), 0.5*zeros(numel(event_ts(event_inds),1)), 'rx')
hold on
plot(ts,new_TTL_file)
ylim([-0.1, 1.1])
