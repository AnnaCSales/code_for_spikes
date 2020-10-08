function [wf_feats]=find_wf_features(path,spikeStruct, iUnit, pl)
%pl is 0 or 1 depending on whether plots are required or not.
% Extract basic info about the waveform from the spikeStruct & plot, if
% needed
fs=spikeStruct.sample_rate;
if pl
    wave_fig=figure('Name', 'Average waveform for each cluster', 'color' ,'w', 'units', 'centimeters', 'pos', [1, 2, 6, 6]);
end
ts_= spikeStruct.timesSorted{iUnit}; 
wave=spikeStruct.av_waveform{iUnit};
errb=spikeStruct.std_waveform{iUnit};
wave_time=1000/fs *(1:length(wave));
c_chan=spikeStruct.c_channel(iUnit); 

if pl
    plot(wave_time, wave)
     hold on
     xlabel('Time, (ms)')
     ylabel('Voltage \muV')
     title([path ', unit: ' num2str(iUnit)]);
end
%%
%The central negative peak is fixed at t~1.4ms. Need to try to get
%information about the width of waveforms either side.

window_width= 1/fs * length(wave); %the width of the entire period, in seconds.
cent_ind=0.5*(length(wave));
cent_peak= 1000/fs * cent_ind; %the time of the central negative peak

%first find the width of the central peak

[cent_amp,cent_ind]=min(wave); %voltage at central peak
if pl
    plot(wave_time(cent_ind),cent_amp, 'ro')
end
%define the baseline voltage - use to define relative heights of peaks
%etc.
bl_v=mean(wave(1:10));
cent_rel_amp=abs(bl_v-cent_amp);

%Value at FWHM, central peak
wid_point=bl_v-(0.5*cent_rel_amp); %FWHM for now, can change to 1/e etc.

right_wid=interp1(wave(cent_ind:cent_ind+30), wave_time(cent_ind:cent_ind+30), wid_point);
left_wid=interp1(wave(cent_ind-10:cent_ind), wave_time(cent_ind-10:cent_ind), wid_point)  ;    
t_width=right_wid-left_wid;  %FWHM width in ms

if pl
    plot([left_wid, right_wid],[wid_point,wid_point], 'x')
end
%%
%find peaks (turning points)
[peak_vals, peak_inds]=findpeaks(wave, 'MinPeakWidth', 2, 'MinPeakDistance', 5,'MinPeakProminence',3) ;

if pl
    plot(wave_time(peak_inds), peak_vals, 'or')
end

%Find the start of the waveform
early_inds=1:0.5*(length(wave))-2;
early=wave(early_inds); %the first half of the waveform             
flat_=mean((diff(early(1:20))));  
flat_sd=std(diff(early(1:20)));

[~,i_st]=find( abs(diff(early))>= abs(flat_)+6*abs(flat_sd));

if length(i_st)==0      
       [~,i_st]=find( abs(diff(early))>= abs(flat_)+5*abs(flat_sd));
end

if pl
    plot(wave_time(early_inds(i_st(1))), wave(early_inds(i_st(1))), 'diamond' )
end

start_wf=wave_time(early_inds(i_st(1)));
%% Deal with ones with missing peaks, so that peaks_ind is always 2-valued.

%1. No first peak:

if ~length(peak_inds(peak_inds<cent_ind))  
    peak_inds=[NaN,peak_inds]; %just stick a NaN there to indicate no peak
end


%now deal with the ones where there's no recognisable peak after the
%  central one - stick in the max. 
if ~length(peak_inds(peak_inds>cent_ind))  %need to look for that first drop in voltage if there's no peak
    late_inds=cent_ind+2:length(wave);
    late=wave(late_inds);    
    [v_max, i_max]=max(late);
    late_max_ind=late_inds(i_max);
    late_max_val=v_max;
    peak_inds(2)=late_max_ind;
    peak_vals(2)=v_max;
    if pl
    plot(wave_time(peak_inds(2)), v_max, 'co');
    end
end

   
%% Now gather together all the info and return it.

wf_feats.wf=wave;
wf_feats.wf_ts=wave_time;
wf_feats.c_chan=c_chan;
wf_feats.cent_ind=cent_ind;
wf_feats.bl_v=bl_v;
wf_feats.FWHM=t_width;
wf_feats.peak_vals=peak_vals;
wf_feats.peak_inds=peak_inds;
wf_feats.start_wf=start_wf;
wf_feats.cent_amp=cent_amp;

if pl
    d_=strsplit(path);
    dt=d_{end};
    saved_name=[dt '_unit_ ' num2str(iUnit)];  %NB PATH HARD CODED!!
     savefig( wave_fig,['E:\New code\fig2\waveform_pics' '\' saved_name '.fig']);  
    close(wave_fig)
end
end