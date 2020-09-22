% This is a script for analysing mouse ECG/respiratory data from the NTS prep. It does the
% following:
%  - reads in openEphys data for mouse ECG,
%  - filters it to highlight the parts of the signal relating to the respiratory cycle / ECG
%  - find the number of heart beats / breaths during a specified period, to
%    calculate bradycardia/pnea relative to markers (e.g. laser
%    activation)
% You need to update the filename fields, e.g. ADC_fn and TTL_fn, so that
% the script can find the ADC data and the TTL events, and also make sure
% that the 'laser_BNC_chan' correctly shows where the laser BNC was plugged
% in to the input/output board.
% Anna Sales, UoB, August 2020.


%% Start by reading in the raw ADC data and the TTL information

ADC_fn='D:\NTS_Pabitra\NTS_new_data\EMG_ECG\100_ADC2.continuous'; %this is a string containing the location of the data file
[raw_sig, ADC_ts_raw,  ADC_info] = load_open_ephys_data(ADC_fn);  %returns the data, the timebase, and a struct with more useful information

% Read in the events file, for the TTL
TTL_fn='D:\NTS_Pabitra\NTS_new_data\EMG_ECG\all_channels.events'
[events, ev_ts_raw, ev_info] = load_open_ephys_data(TTL_fn);

% Set timestamps so that zero is at the start of the recording (as for the
% spiking data from kilosort/PHY)
ADC_tbase=ADC_ts_raw-ADC_ts_raw(1); % timebase should start at zero
fs=ADC_info.header.sampleRate;  %read the sampling rate
ev_ts=ev_ts_raw-ADC_ts_raw(1);  %set to a zero'd starting point
laser_BNC_chan=4; %where the laser TTL BNC was plugged in to the in/out board
evind=find(events==4);  % pull out the laser TTLs
laser_ts_all=ev_ts(evind);

%% Detect large gaps in the rec

blocks=block_detector(ADC_tbase);

%% optional - specify a time period (otherwise use entire recording - comment out as req)

analysis_window=[ADC_tbase(1), 1050]  %the time window we want to use

% analysis_window=[230, 240]; %specify the period we want to analyse, in seconds

% now cut the data, so that we only work with data in the time range
% specified in 'analysis window'

% the indices of data in the specified range:
keep_=find(ADC_tbase>analysis_window(1) & ADC_tbase<analysis_window(2));

%the timebase and raw signal in the time window:
tbase_win=ADC_tbase(keep_);
raw_sig=raw_sig(keep_);

laser_inc=find(laser_ts_all>analysis_window(1) & laser_ts_all<analysis_window(2));
laser_ts=laser_ts_all(laser_inc);


%% Filter, find peaks representing heart beat/ breathing

[ecg, resps]=return_filtered_data(raw_sig, fs);

% find peaks in the data and mark on the plot (may need to play with the
% parameters in the function call)
[pks_resp,times_resp] = findpeaks(resps,fs, 'MinPeakProminence',0.08,'MinPeakDistance',0.1);
times_resp=times_resp+analysis_window(1); %set times relative to start of window

% store a continuous log of gaps between breaths
resp_gaps=[0; diff(times_resp)];
times_gaps=[times_resp-[0;0.5*diff(times_resp)]];  %time points for these gaps (halfway between breaths)

% find the negative going peaks in the ECG signal
[pks_ecg,times_ecg] = findpeaks(-1*ecg,fs, 'MinPeakProminence',0.1,'MinPeakDistance',0.05);
times_ecg=times_ecg+analysis_window(1); %set times relative to start of window

% store the continuous rate (per min) as 60/(time between consecutive beats)
ecg_rate=[0;60./diff(times_ecg)];

%% Now plot 20 secs of data (helps to sanity check, and
% fine tune parameters for peak finding above)

filt_fig=figure('Color', 'w')
subplot(2,1,1)
in_win=find(tbase_win>analysis_window(1) & tbase_win<analysis_window(1)+10);
plot(tbase_win(in_win),resps(in_win), 'b');
ylim([-0.2,0.3]);
ylabel('Filtered signal, \muV');
xlim([analysis_window(1), analysis_window(1)+10]);
xlabel('Time (s)');
hold on
plot([laser_ts, laser_ts], [-0.2, 0.3], 'Color', [1 0 0 0.1])
plot(times_resp, pks_resp, 'kX')
title('Example data: respiratory signal', 'FontWeight', 'normal')


subplot(2,1,2)
plot(tbase_win(in_win),ecg(in_win), 'b');
ylim([-0.2,0.2]);
ylabel('Filtered signal, \muV');
title('ECG signal','FontWeight', 'normal')
hold on
plot(times_ecg, -1*pks_ecg, 'kX')
xlim([analysis_window(1), analysis_window(1)+10]);
plot([laser_ts, laser_ts], [-0.2, 0.3], 'Color', [1 0 0 0.1])
%% Go through the data counting peaks  - get bradycardia/npea in specified periods.

% First pull out the start and end of laser blocks (not the start/stop of
% individual pulses).  We assume these are event separated by >20s
% The widget processes the individual TTLs into groups,
% calculating freq, pulse duration, number of pulses
% stims: table with cols  time of event, freq.
% stims_by_type: same, but grouped by event type.'
% utils: labels and markers for plotting

[stims, stims_by_type, stim_utils]=laserTTLwidget(laser_ts, 20);

%the following block makes some markers
for h=1:size(stims_by_type,1)
    stim_=stims_by_type{h}
    %extract info about this stimtype
    freq=stim_(1, 2)
    num=stim_(1, 4)
    len=stim_(1, 3)
    
    %create some markers for use in plotting
    if freq==0
        stimMarkers{h}=[0; len/1000];
    else
        starts=(0:num-1)* 1/freq;
        ends=starts+len/1000;  %length is in ms in the table
        stamps=zeros(2*length(starts), 1);
        stamps(1:2:end)=starts;
        stamps(2:2:end)=ends;
        stimMarkers{h}=stamps;
    end
    
    %make some labels to describe each type of laser event
    if freq==0
        str1='Single pulse';
        stim_type_labels{h}=[str1 ', ' num2str(stims_by_type(h,2)) ' ms'];
    else
        stim_type_labels{h}=[num2str(num) ' pulses, ' num2str(len) ' ms, '  num2str(freq) ' hz'];
    end
    
    
end
%%
%period to use for calculating mean rates before and after
comp_period=10; %start at 10s

for types=1:length(stims_by_type)
    this_type=stims_by_type{types}  %all events of this type
    %     this_type=this_type{1};
    
    num_=this_type(1, 4)
    freq_=this_type(1, 2)
    pulse_dur=this_type(1,3)
    event_duration=1/freq_ * num_ ;  % the total time for one of these events
    
    resp_results=[]; % for storing results.
    ecg_results=[];
    
    for stim_num=1:length(this_type(:,1)) %now go through the individual events
        
        start_=this_type(stim_num,1);
        end_=start_ + event_duration;
        
        % find the times of peaks in the following periods:
        % baseline is 2 comp_periods back from start
        % pre_stim is 1 comp_period back from start
        % post_stim is 1 com_period after end
        % recovery is 2 comp_perods after end
        % periods={'baseline', 'pre', 'opto','post', 'recovery' }
        cp_pos=[2,1,0,1,2]; %positions of comparator periods.
        for period=1:5;
            t=cp_pos(period);
            % take all the event during a given period, plus the two either
            % side (so that we can get an instantaneous rate for all
            % events)
            if t==0  % if it's during the stim
                % indices of ecg/resp peaks during the stim period:
                resp_inds= find(times_resp>start_ & times_resp <= end_) ;
                ecg_inds= find(times_ecg>start_ & times_ecg <= end_) ;
            elseif period<3
                resp_inds = find(times_resp>=(start_-  t*comp_period) & times_resp < (start_- (t-1)*comp_period) ) ;
                ecg_inds = find( times_ecg>=(start_-  t*comp_period) & times_ecg < (start_- (t-1)*comp_period) ) ;
            elseif period>3
                resp_inds= find( times_resp> (end_+ (t-1)*comp_period) & times_resp < (end_+ t*comp_period) ) ;
                ecg_inds= find( times_ecg> (end_+ (t-1)*comp_period) & times_ecg < (end_+ t*comp_period) ) ;
            end
            
            %           resp_times=[times_resp(resp_inds(1)-1); times_resp(resp_inds);times_resp(resp_inds(end)+1) ];
            %           ecg_times=[times_ecg(ecg_inds(1)-1); times_ecg(ecg_inds);times_ecg(ecg_inds(end)+1) ]  ;
            
            resp_times=[ times_resp(resp_inds) ];
            ecg_times=[ times_ecg(ecg_inds) ]  ;
            
            mean_resp=mean(diff(resp_times));  %work with time intervals between breaths.
            min_resp=min(diff(resp_times));
            max_resp=max(diff(resp_times));
            mean_ecg=mean(60./(diff(ecg_times)));  %work with heart rate.
            min_ecg=min(60./(diff(ecg_times)));
            max_ecg=max(60./(diff(ecg_times)));
            
            resp_results(period, :, stim_num)=[mean_resp, min_resp, max_resp];
            ecg_results(period,:, stim_num)=[mean_ecg, min_ecg, max_ecg];
            
        end
        
    end
    all_resp_results{types}=resp_results;
    all_ecg_results{types}=ecg_results;
end
%% Write the raw data to xls form
cell_labels_c={'Mean (bpm)', 'Min (bpm)', 'Max (bpm)'};
time_labels={'baseline', 'pre', 'opto', 'post', 'recovery'}';
fn_c='bradycardia.xls';
desc_c={'Heart rate:'};
for t=1:length(all_ecg_results)
    this_set=all_ecg_results{t};
    stim_label={stim_type_labels{t}};
    
    for s=1:size(this_set,3)
       this_sheet=this_set(:,:,s);
       xlswrite(fn_c,cell_labels_c,s,'B1:D1');
       xlswrite(fn_c, time_labels, s,'A2:A6');
       xlswrite(fn_c, this_sheet, s, 'B2:D6');
       xlswrite(fn_c, stim_label, s, 'F1');
       xlswrite(fn_c, desc_c, s, 'A1');
    end
    
end


fn_r='bradypnoea.xls';
cell_labels_r={'Mean (s)', 'Min (s)', 'Max (s)'};
desc_r={'Interbreath interval:'};
for t=1:length(all_resp_results)
    this_set=all_resp_results{t};
    stim_label={stim_type_labels{t}};
    for s=1:size(this_set,3)
       this_sheet=this_set(:,:,s);
       xlswrite(fn_r,cell_labels_r,s,'B1:D1');
       xlswrite(fn_r, time_labels, s,'A2:A6');
       xlswrite(fn_r, this_sheet, s, 'B2:D6');
       xlswrite(fn_r, stim_label, s, 'F1');
       xlswrite(fn_r, desc_r, s, 'A1');
    end
    
end
%% Now calculate stats as per powerpoint.
%Tables arranged mean, min, max across and event number down
type_=1;

stim_paras=(stims_by_type{type_})
num_=stim_paras(1, 4)
freq_=stim_paras(1, 2)
pulse_dur=stim_paras(1,3)
event_duration=1/freq_ * num_ ;


laserFig      = figure('color','w','NumberTitle','off', 'name',stim_type_labels{type_}, 'units', 'centimeters', 'pos',[5 2 24 17]);
laserTabGroup = uitabgroup(laserFig,'TabLocation','Left');

ecg_res=all_ecg_results{1};
resp_res=all_resp_results{1};
brady_c=[];
brady_p=[];
for t=1:size(resp_res,3) %number of events
    % calculate key stats for each event
    event_res=resp_res(:,:,t);
    event_ecg=ecg_res(:,:,t);
    
    baseline_bradyc=((event_ecg(1,1)-event_ecg(2,2))/event_ecg(1,1)) * 100;
    opto_bradyc=((event_ecg(1,1)-event_ecg(3,2))/event_ecg(1,1)) * 100;
    recovery_bradyc= ((event_ecg(1,1)-event_ecg(5,2))/event_ecg(1,1)) * 100;
    delta_bradyc= opto_bradyc- baseline_bradyc;
    
    baseline_bradyp=((event_res(2,3)-event_res(1,1))/event_res(1,1)) * 100;
    opto_bradyp=((event_res(3,3)-event_res(1,1))/event_res(1,1)) * 100;
    recovery_bradyp=((event_res(5,3)-event_res(1,1))/event_res(1,1)) * 100;
    delta_bradyp= opto_bradyp- baseline_bradyp;
    
    brady_c(t,:)=[baseline_bradyc, opto_bradyc, recovery_bradyc, delta_bradyc]
    brady_p(t,:)=[baseline_bradyp, opto_bradyp, recovery_bradyp, delta_bradyp]
    
    
    %plot each event as a sanity check
    figure(laserFig)
    laser_tab_sp = uitab(laserTabGroup, 'Title', ['Stim #: ' num2str(t) ],'BackgroundColor','w');
    axes('Parent',laser_tab_sp);
    
    keep_=find(tbase_win>stim_paras(t,1)-20 & tbase_win<stim_paras(t,1)+20);
    
    sb_resp1=subplot('Position',[0.13, 0.74, 0.7750, 0.23]);
    sb_resp2=subplot('Position', [0.13, 0.58, 0.7750, 0.11]);
    sb_ecg1=subplot('Position', [0.13, 0.23, 0.7750, 0.23]);
    sb_ecg2=subplot('Position', [0.13, 0.06, 0.7750, 0.11]);
    
    subplot(sb_resp1)
    plot(tbase_win(keep_),resps(keep_), 'b');
    ylim([-0.1,0.2]);
    hold on
    plot(times_resp, pks_resp, 'ko')
    xlim([min(tbase_win(keep_)), max(tbase_win(keep_))])
    ylabel('Filtered signal, \muV');
    hold on
    x_patch= [stim_paras(t,1), stim_paras(t,1)+event_duration, stim_paras(t,1)+event_duration, stim_paras(t,1)];
    y_patch=[-0.2, -0.2, 0.3, 0.3];
    patch(x_patch, y_patch, 'r','FaceAlpha',.4, 'EdgeColor', 'none');
    title('Respiratory signal', 'Fontweight', 'normal')
    xlim([min(tbase_win(keep_)), max(tbase_win(keep_))])
    box off
    
    subplot(sb_resp2)
    keepR=find(times_gaps>stim_paras(t,1)-20 & times_gaps<stim_paras(t,1)+20);
    plot(times_gaps(keepR), resp_gaps(keepR), 'Color', rgb('MediumSeaGreen'))
    ylabel({'Interbreath' 'interval (s)'})
    xlabel('Time (s)');
    xlim([min(tbase_win(keep_)), max(tbase_win(keep_))])
    box off
    hold on
    y_patch=[-5, -5, 5, 5];
    patch(x_patch, y_patch, 'r','FaceAlpha',.4, 'EdgeColor', 'none');
    ylim([0.5,2.5])
    text(max(tbase_win(keep_))-10,2.4, ['\Delta bradypnoea during opto = ' num2str(brady_p(t,4)), '%'], 'FontSize', 8 )
    
    subplot(sb_ecg1)
    plot(tbase_win(keep_),ecg(keep_), 'b');
    hold on
    xlim([min(tbase_win(keep_)), max(tbase_win(keep_))])
    ylabel('Filtered signal, \muV');
    hold on
    x_patch= [stim_paras(t,1), stim_paras(t,1)+event_duration, stim_paras(t,1)+event_duration, stim_paras(t,1)];
    y_patch=[-0.5, -0.5, 0.5, 0.5];
    patch(x_patch, y_patch, 'r','FaceAlpha',.2, 'EdgeColor', 'none');
    title('ECG signal', 'Fontweight' , 'normal')
    ylim([-0.1,0.2]);
    xlim([min(tbase_win(keep_)), max(tbase_win(keep_))])
    box off
    
    subplot(sb_ecg2)
    keepE=find(times_ecg>stim_paras(t,1)-20 & times_ecg<stim_paras(t,1)+20);
    plot(times_ecg(keepE), ecg_rate(keepE), 'Color', rgb('MediumSeaGreen'))
    ylabel({'Heart rate' 'per min'})
    box off
    hold on
    aa=gca;
    ylim([370, 420])
    y_patch=[aa.YLim(1),aa.YLim(1),aa.YLim(2),aa.YLim(2)];
    patch(x_patch, y_patch, 'r','FaceAlpha',.4, 'EdgeColor', 'none');
    xlim([min(tbase_win(keep_)), max(tbase_win(keep_))])
    xlabel('Time (s)')
    ylim([370, 420])
    text(max(tbase_win(keep_))-10,aa.YLim(2), ['\Delta bradycardia during opto = ' num2str(brady_c(t,4)) '%'], 'FontSize', 8 )
    
end

%mean delta_opto for all events of this type
mean_brady_c=mean(brady_c(:,4), 1)
sem_brady_c=nansem(brady_c(:,4),1)
% [h,p]=ttest(brady_c(:,4));
[p_c,~,~]=signrank(brady_c(:,4))

mean_brady_p=mean(brady_p(:,4), 1)
sem_brady_p=nansem(brady_p(:,4),1)
% [h,p]=ttest(brady_p(:,4))
[p_p,~,~]=signrank(brady_p(:,4))


