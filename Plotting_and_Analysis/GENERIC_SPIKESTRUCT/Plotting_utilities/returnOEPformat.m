function [file_prefix, file_postfix] = returnOEPformat(datapath)
% This function takes a directory name, examines the open Ephys files in
% that directory and returns the filename format. This is useful for
% specifying a directory name and automatically reading in all the data in
% that folder.

file_list=dir(datapath);
chk4cont = regexp({file_list.name}, '.continuous', 'once');
cont_matches=find(~cellfun(@isempty, chk4cont));

if length(cont_matches)>0
    
    fprintf('\n Continuous data...reading example file....\n')
    ex_file=file_list(cont_matches(1)).name;
    %only read the first sample - this will give us the start time,
    %and metadata
    [ex_data, ts_start, info] = load_open_ephys_data([datapath, ex_file],'Indices', 1);   
    name_parts=strsplit(ex_file, '_');
    file_prefix=[name_parts{1} '_'];%the bits around 'CH1', i.e. 113_CH1_2, pre is 113_, post is _2
    
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
end

end

