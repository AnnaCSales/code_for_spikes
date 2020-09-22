function [stims, stims_by_type, stim_utils] = laserTTLwidget(laser_ts, varargin)
% If the TTLs on a given channel represent the start and end of laser
% flashes, this code parses the TTLs into frequency, pulse duration, and
% number of flashes. Takes the TTL channel on which the laser is represented, and an
% optional argument for the rough separation between laser events (i.e. 
% for an argument 't', anything  more than t seconds apart is designiated 
% as a separate event.

% Returns:

% 'stims', a matrix with column headings
% [time of first flash, frequency, pulse duration, number of flashes]

% 'stims_by_type', a cell array where each cell holds the times of
% stimulations of identical type, in the same format as above.

if ~isempty(varargin)
    gaps=varargin{1};  %gaps is the minimum time between laser blocks.
else
    gaps=2
end

laser_ts=unique(laser_ts);  
lts_gaps=diff(laser_ts);  %index n of lts_gaps is the gap between TTL n+1 and n
[l_inds, ~]=find(lts_gaps>gaps); %pulls out the gaps between laser TTLs which are more than 2s apart
stim_starts=[1; (l_inds+1)];
num_stims=length(stim_starts);
stims=zeros(length(stim_starts), 4) ; %this will hold details of all of the pulse trains in this recording.

%extract all of the laser pulse trains and store information about them
for p=1:num_stims
    
    st_ind=stim_starts(p);
    
    if p==length(stim_starts)
        end_ind=(length(laser_ts));
    else
        end_ind=stim_starts(p+1)-1;
    end
    
    num_=(end_ind-st_ind+1) /2;
    length_=laser_ts(st_ind+1)-laser_ts(st_ind);
    if (end_ind-st_ind==1)
        freq_=0;  %for single pulse stimulations
    else
        freq_= round( (1 / ( laser_ts(st_ind+2)-laser_ts(st_ind))) , 3); %freqs.
    end
    
    num_=round(num_, 1);
    length_=round(length_, 3)* 1000 ;%in ms
    
    stims(p,:)=[laser_ts(st_ind), freq_, length_, num_];
    %table with times relative to zero'd timestamps, then freq, length, num
end


laserTTL=stims(:,1);

%do a bit more work with the laser timestamps to pull out indicies of
%identical events
[stimtypes, ~,r_inds]=unique(stims(:,[2:4]), 'rows', 'first');

stims_by_type={};
 for h=1:size(stimtypes,1)
     
     stims_by_type{h}=stims(find(r_inds==h), :); %find rows of that type, store
     
     %extract info about this stimtype
     freq=stimtypes(h, 1);
     num=stimtypes(h, 3);
     len=stimtypes(h, 2);
     
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
     if stimtypes(h,1)==0
         str1='Single pulse';
         stim_type_labels{h}=[str1 ', ' num2str(stimtypes(h,2)) ' ms'];
     else
         stim_type_labels{h}=[num2str(stimtypes(h,3)) ' pulses, ' num2str(stimtypes(h,2)) ' ms '  num2str(stimtypes(h,1)) ' hz'];
     end
        
 end
stim_utils.Markers=stimMarkers;
stim_utils.labels=stim_type_labels;
 
end

