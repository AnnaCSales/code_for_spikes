function [TTLs] = returnTTLs(dirname,varargin)
% Code takes either continuous or binary OEP data and returns TTL input on 
% each of the 8 digital channels, alongside the manually entered information.
% Expects input arguments (strings) in the form 
% TTLs= returnTTLs('C:\exampledir') ; for OEP continuous recordings
% NB there must be an example .continuous file in the folder to correctly
% extract the TTLs!
% Varargin provides options for binary data (still being developed)


% Anna Sales 31/07/2020

% Work out what files are in the directory given
file_list=dir(dirname);
chk4cont = regexp({file_list.name}, 'continuous', 'once');
cont_matches=find(~cellfun(@isempty, chk4cont));

if length(cont_matches)>0
    
    fprintf('\n Continuous data...reading example file....\n')
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


    [data_raw, ts_raw,  info] = load_open_ephys_data([dirname ex_file]);
    [event_data, event_ts, event_info] = load_open_ephys_data([dirname ['all_channels' file_postfix '.events']]);
    manual_events=fileread([dirname ['messages' file_postfix '.events']]);
    fs=info.header.sampleRate;
    
    %parse the events file:
    find_newlines=regexp(manual_events, '\n');
    messages={};
    if ~isempty(find_newlines)
        for n=2:length(find_newlines)  %ignore the first one as it's always an internal messages
            if n<length(find_newlines)
              messages{n-1}=manual_events(find_newlines(n):find_newlines(n+1)-1);
            else
              messages{n-1}=manual_events(find_newlines(n):end);
            end
        end
        messages(end)=[]; %last one is nonsense also.
    end
    %Check the recording isn't split into blocks.

    ts_gaps=diff(ts_raw);  %index n of lts_gaps is the gap between TTL n+1 and n
    [l_inds, ~]=find(ts_gaps>1); %pulls out the gaps between timestamps that are are more than 1s apart
    block_times(:,1)=ts_raw([1; (l_inds+1)]);  %first col is starts
    block_times(:,2)=ts_raw( [l_inds; length(ts_raw)]);   %second col is ends
    num_block=length(block_times(:, 1));
    
    if num_block>1
        fprintf('This recording is in multiple blocks -  exiting');
    else
        t_first=ts_raw(1); %time of first timestamp - need to deduct this from all TTLs to get times ref to zero.

        %  Extract TTLs 1-8
        TTLtimes=cell(1,8); 
        for iTTL=1:8       
            TTL_ts=unique(event_ts( find(event_data==iTTL-1)));  %find of elements which have footshock timestamp.
            TTL_ts=TTL_ts-t_first;  %offset to start at zero, like the PHY output.
            TTLtimes{iTTL}=TTL_ts;
        end

        % now do the manual events
        TTL_times=[];
        TTL_labels={};

        for m=1:length(messages)
                [time, label]=strtok(messages{m}, ' ');
                if ~isempty(label)
                  TTL_times(m)=(str2num(time)* (1/fs));
                  TTL_labels{m}=label;
                end
        end

            manTTL.TTL_times=TTL_times-t_first;
            manTTL.TTL_labels=TTL_labels;
            
            
            TTLs.digital=TTLtimes;
            TTLs.manual=manTTL;
            save([dirname 'TTLs'], 'TTLs');
    end    
else
    
    fprintf('\n No continuous files found, looking for binary and npy files.')  
    try
%         fullpath=fullfile(dirname);
%         pathparts=strsplit(fullpath, '\');
%         cont_ind= (regexp(pathparts, 'continuous', 'once')); %find the folder marked 'continuous'
%         matches=find(~cellfun(@isempty, cont_ind));
%         path_to_structurefile=[strjoin(pathparts(1:matches-1),'\' ) '\'];
        
        path_to_structurefile=varargin{1};
        events_data=load_open_ephys_binary([path_to_structurefile 'structure.oebin'], 'events',1);
        fs=events_data.Header.sample_rate;   
        data_ts=readNPY([dirname 'timestamps.npy']); %timestamps for the voltage data.
        
        
        % get digital TTLs
        TTL_ts_=events_data.Timestamps; %times of all TTLs
        TTL_chans=events_data.ChannelIndex;  %chan corresponding to each timestamp.
  
        %reference all to zero
        TTL_ts=double(TTL_ts_-data_ts(1))/fs;

        TTL_times={}; %pull out individual channels
        for iTTL=1:8
            thisTTL= find(TTL_chans==iTTL);
            TTLtimes{iTTL}=TTL_ts(thisTTL); %timestamps for TTLs on this channel.
        end

        %now do the manual ones.
        
        json=jsondecode(fileread([path_to_structurefile 'structure.oebin']));
        
        [rhythm, message_centre]=json.events.folder_name; %paths for internal/external messages
        path_to_events=[path_to_structurefile 'events\' message_centre]
        messages = readNPYtext([path_to_events 'text.npy']);
        message_ts=readNPY([path_to_events 'timestamps.npy']);

        manTTL.TTL_times=double(message_ts)/fs;
        manTTL.TTL_labels=messages;

        TTLs.digital=TTLtimes;
        TTLs.manual=manTTL;
        save([dirname 'TTLs'], 'TTLs');
    catch
        fprintf('\nCannot process TTL data in either binary or continuous format, exiting\n')
        fprintf('Check folder names, directory names and contents\n')
    end
end

