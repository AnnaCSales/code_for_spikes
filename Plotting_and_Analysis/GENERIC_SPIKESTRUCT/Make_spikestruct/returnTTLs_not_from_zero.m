function [TTLs] = returnTTLs_not_from_zero(dirname)
% This function takes continuous OEP data and returns TTL input on 
% each of the 8 digital channels, alongside the manually entered information.
% Timestamps are NOT relative to a zero starting point so should not be
% used with clustered data from PHY.
% Expects input arguments (strings) in the form. 
% 'dirname' is the path to the folder containing the OEP data.

% Example:
% TTLs= returnTTLs('C:\exampledir') ; for OEP continuous recordings
% NB there must be an example .continuous file in the folder to correctly
% extract the TTLs.

% Anna Sales 31/07/2020


% Work out what files are in the directory given, and the format of their
% names
file_list=dir(dirname);
chk4cont = regexp({file_list.name}, '.continuous', 'once');
cont_matches=find(~cellfun(@isempty, chk4cont));

if length(cont_matches)>0
    
    fprintf('\n Continuous data found...reading example file....\n')
    ex_file=file_list(cont_matches(1)).name;
    name_parts=strsplit(ex_file, '_');
    prefix=[name_parts{1} '_'];%the bits around 'CH1', i.e. 113_CH1_2, pre is 113_, post is _2
    
    if numel(name_parts)>2
        post_=strsplit(name_parts{3}, '.');
        
        if numel(post_)>1
            file_postfix=['_' post_{1}];
        else
            file_postfix=[];
        end
    else
        file_postfix=[];
    end
    
    % Load in some example data and the events file
    [data_raw, ts_raw,  info] = load_open_ephys_data([dirname ex_file]);
    [event_data, event_ts, event_info] = load_open_ephys_data([dirname ['all_channels' file_postfix '.events']]);
    manual_events=fileread([dirname ['messages' file_postfix '.events']]);
    fs=info.header.sampleRate;
    
    %parse the events file (this is a text file containing manually entered info):
    find_newlines=regexp(manual_events, '\n');
    messages={};
    for n=2:length(find_newlines)  %ignore the first one as it's always an internal messages
        if n<length(find_newlines)
            messages{n-1}=manual_events(find_newlines(n):find_newlines(n+1)-1);
        else
            messages{n-1}=manual_events(find_newlines(n):end);
        end
    end
    messages(end)=[]; %last one is nonsense also.
    
    %Check the recording isn't split into blocks.
    ts_gaps=diff(ts_raw);  %index n of lts_gaps is the gap between TTL n+1 and n
    [l_inds, ~]=find(ts_gaps>1); %pulls out the gaps between timestamps that are are more than 1s apart
    block_times(:,1)=ts_raw([1; (l_inds+1)]);  %first col is starts
    block_times(:,2)=ts_raw( [l_inds; length(ts_raw)]);   %second col is ends
    num_block=length(block_times(:, 1));
    
    if num_block>1
        fprintf('WARNING -  this recording is in multiple blocks ');
    else
        
        %  Extract TTLs 1-8
     TTLtimes=cell(1,8);
     for iTTL=1:8
         TTL_ts=unique(event_ts( find(event_data==iTTL-1)));  %find of elements which have footshock timestamp.
         TTL_ts=TTL_ts;  %offset to start at zero, like the PHY output.
         TTLtimes{iTTL}=TTL_ts;
     end
        
        % now do the manual events
     TTL_times=[];
     TTL_labels={};
        
     for m=1:length(messages)
         [time, label]=strtok(messages{m}, ' ');
         TTL_times(m)=(str2num(time)* (1/fs));
         TTL_labels{m}=label;
     end
        
     manTTL.TTL_times=TTL_times;
     manTTL.TTL_labels=TTL_labels;
        
        
     TTLs.digital=TTLtimes;
     TTLs.manual=manTTL;
        save([dirname 'TTLs'], 'TTLs');
    end
else    
    fprintf('\n No continuous files found.')  
end

