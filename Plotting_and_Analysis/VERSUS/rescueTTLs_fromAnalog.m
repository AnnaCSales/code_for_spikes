[event_data_all, event_ts_all, event_info] = load_open_ephys_data('C:\Versus data\281020\2020-10-28_16-13-14\all_channels.events');
[data_analog, ts_, info] = load_open_ephys_data('C:\Versus data\281020\2020-10-28_16-13-14\100_ADC7.continuous');   
[data_ch_, ~, info] = load_open_ephys_data('C:\Versus data\281020\2020-10-28_16-13-14\CAR\102_CH19.continuous');   

%% Try to extract TTLs from artefact on the raw data
% Use the initial positive deflection above 500uv.
data_ch=abs(data_ch_);
high_inds=find(data_ch>350);
%Look for blocks above 350 separated by at least 0.4s (fastest we went was
%2Hz)
gaps=find(diff(high_inds)> (fs*0.4));
figure
plot(ts_, data_ch_)
hold on
% plot(ts_(high_inds), data_ch_(high_inds), 'r')
plot(ts_(high_inds(gaps+1)), 500, 'go')

footshock_TTLs=ts_(high_inds(gaps+1));
footshock_TTLs_zero=footshock_TTLs-ts_(1);