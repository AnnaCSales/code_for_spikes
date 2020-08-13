function [stims, stims_by_type] = laserTTLwidget(laser_ts)
% If the TTLs on a given channel represent the start and end of laser
% flashes, this code parses the TTLs into frequency, pulse duration, and
% number of flashes (anything >2 seconds apart is designiated as a separate
% stim). Takes the TTL channel on which the laser is represented, returns:

% 'stims', a matrix with column headings
% [time of first flash, frequency, pulse duration, number of flashes]

% 'stims_by_type', a cell array where each cell holds the times of
% stimulations of identical type, in the same format as above.


lts_gaps=diff(laser_ts);  %index n of lts_gaps is the gap between TTL n+1 and n
[l_inds, ~]=find(lts_gaps>2); %pulls out the gaps between laser TTLs which are more than 2s apart
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
        freq_= round( (1 / ( laser_ts(st_ind+2)-laser_ts(st_ind))) , 0); %integer freqs.
    end
    
    num_=round(num_, 1);
    length_=round(length_, 2)* 1000 ;%in ms
    
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
           
 end

 
end

