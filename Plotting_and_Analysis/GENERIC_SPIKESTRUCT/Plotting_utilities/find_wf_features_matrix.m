function [amps, FWHMs]=find_wf_features_matrix(spike_data, fs, pl)
%Extract basic features of spikes, using a matrix of spike waveforms
%(spike_data, where each row is the trace for an individual spike. 
% fs is the sampling rate.
% pl is an optional value to plot five example waveforms with their
% features.

n_spikes=size(spike_data, 1);
amps=zeros(n_spikes, 1);
FWHMs=zeros(n_spikes,1);
if pl
    shuffled_spikes=randperm(n_spikes);
    spk_examples=shuffled_spikes(1:5); % five randomly chosen spikes to plot.
    wave_fig=figure('color' ,'w', 'units', 'normalized', 'pos', [0.1, 0.1, 0.7, 0.3]);
end
rw=[];
n_wave_samples=size(spike_data,2); % number of samples per wave
wave_time=1000/fs *(1:n_wave_samples);

%The central negative peak is fixed at t~1.4ms. Need to try to get
%information about the width of waveforms either side.

window_width= 1/fs * n_wave_samples; %the width of the entire period, in seconds.
cent_ind=0.5*n_wave_samples;
cent_peak= 1000/fs * cent_ind; %the time of the central negative peak

p=1;

for s=1:size(spike_data,1)
    
    try
        plotflag=0;
        wave=spike_data(s,:);
      

        if ismember(s, spk_examples) & pl==1
            figure(wave_fig)
            subplot(1,5,p)
            plot(wave_time, wave)
            title(['Example #' num2str(p)], 'FontWeight', 'normal', 'FontSize', 8)
            hold on
            xlabel('Time, (ms)')
            ylabel('Voltage \muV')
            xlim([0, wave_time(end)])
            p=p+1;
            plotflag=1;
        end

        %%
        %The central negative peak is fixed at t~1.4ms. Need to try to get
        %information about the width of waveforms either side.

        %first find the width of the central peak

        [cent_amp,cent_ind]=min(wave); %voltage at central peak
        if plotflag
            plot(wave_time(cent_ind),cent_amp, 'ro')
        end
        % %define the baseline voltage - use to define relative heights of peaks
        % %etc.
        bl_v=mean(wave(1:30));
        cent_rel_amp=abs(bl_v-cent_amp);

        %Value at FWHM, central peak
        wid_point=bl_v-(0.5*cent_rel_amp); %FWHM for now, can change to 1/e etc.
        rise_num=5;
        try
            right_wid=interp1(wave(cent_ind:cent_ind+rise_num), wave_time(cent_ind:cent_ind+rise_num), wid_point); 
            left_wid=interp1(wave(cent_ind-rise_num:cent_ind), wave_time(cent_ind-rise_num:cent_ind), wid_point)  ;
        catch
            rise_num=2;
            right_wid=interp1(wave(cent_ind+1:cent_ind+rise_num), wave_time(cent_ind+1:cent_ind+rise_num), wid_point); 
            left_wid=interp1(wave(cent_ind-rise_num:cent_ind-1), wave_time(cent_ind-rise_num:cent_ind-1), wid_point)  ;
        end
            
            t_width=right_wid-left_wid;  %FWHM width in ms
            rw(s)=right_wid;
        if plotflag
            plot([left_wid, right_wid],[wid_point,wid_point], 'x')
        end

        %% Now gather together all the info and return it.

        amps(s)=-cent_rel_amp;
        FWHMs(s)=t_width;
    catch
      amps(s)=NaN;
      FWHMs(s)=NaN;
      
    end   
end
fprintf('\nCould not process: %d waveforms \n', sum(isnan(amps)))
