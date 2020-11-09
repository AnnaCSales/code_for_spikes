function info = make_bin_in_blocks()
% Import and chunking of Open Ephys data into equal length .bin
% files, followed by merging into one overall binary file. Creates files
% suitable for clustering using Kilosort, from separate .continuous files.
% Navigate to the folder containing the .continuous files to be merged, and
% then run the function. 
% NB this is currently set up for 32 channels. The variable 'nChans' will
% need to be updated if a different number is needed.

%path to OEP utilities folders:
% addpath(genpath('D:\Code\MATLAB\Packages\OEPanalysis-tools-master'))
% addpath(genpath('D:\Code\MATLAB\Packages\npy-matlab-master'))

PathName=pwd;%uigetdir;

input_directory=[PathName];
output_directory = input_directory;
nChans=32;

%pull out the naming format of the files in the specified directory.
directory_contents = dir(input_directory);
contPattern = '*.continuous';
isSettingsFile = regexp({directory_contents.name}, 'continuous', 'once');
matches=find(~cellfun(@isempty, isSettingsFile));
test_name=directory_contents(matches(1)).name;
test_name_parts=strsplit(test_name, '_');

file_prefix=[test_name_parts{1} '_'];%the bits around 'CH1', i.e. 113_CH1_2, pre is 113_, post is _2

if numel(test_name_parts)>2
    post_=strsplit(test_name_parts{3}, '.');
    if numel(post_)>1
        file_postfix=['_' post_{1}];
    else
        file_postfix=[];
    end
else
    file_postfix=[];
end

testfn=[file_prefix 'CH1' file_postfix '.continuous'];

% Specify the bin file length in minutes here:
BlkLengthMin  = 10; 
info = get_session_infoAnna(input_directory, file_postfix);

%% Brute force channel import 
%N.B. 10 minutes of 32 channel recording @ 16bit/30kHz  ~ 1.15GB

ChSaved = 1:nChans;
processor_index = 0;

% grab a channel to get total recording length
fname_in{1} = [input_directory filesep testfn];  %CAN AMEND THIS FOR RECORDINGS WITH _2 etc
evalc('[~,~, fileInfo] = load_open_ephys_data(fname_in{1});'); % Supress commandline output w/ evalc

fLength       = fileInfo.header.bufferSize * length (fileInfo.ts);   % File length in samples
BlkLength     = BlkLengthMin * 60 * fileInfo.header.sampleRate;      % Block length in samples
BlksizeGB     = (BlkLength * 16*  numel(ChSaved))/(1e9*8);           % Block size in GB
noBlks        = ceil(fLength/BlkLength);                             % Number of data blocks to process

fprintf('>>> Recording is %d seconds long (%.2f hours).\n',round(fLength/fileInfo.header.sampleRate),round(fLength/fileInfo.header.sampleRate)/3600 )
fprintf('>>> Writing data as %d blocks of %d minutes, each block is %.3f GB.\n',noBlks,BlkLengthMin,BlksizeGB)

% Cache filenames
for iCh =  1:length(ChSaved)
    fname_in{iCh} = [input_directory filesep [file_prefix 'CH'] int2str(ChSaved(iCh)) file_postfix '.continuous'];
end
for iBlk = 1:noBlks % deal with full blocks first
    fname_out{iBlk} = [get_full_path(output_directory) filesep 'Raw_Part_' num2str(iBlk) '.bin'];
end

% Loop across time in blocks
for iBlk = 1:noBlks
    t1 = tic;
    
    fprintf('>>> Starting Block %d of %d.\n', iBlk, ceil(fLength/BlkLength))
    
    % Preallocate block data (special case for last block as will be time remainder)
    if iBlk ~= noBlks
        for iCh = 1 : length(ChSaved)
            data{iCh} = zeros(BlkLength,1);
            
        end
        
    else
        for iCh = 1 : length(ChSaved)
            data{iCh} =  zeros(fLength-BlkLength*(noBlks-1),1);
        end
        
    end
    flag_ = cell(1,length(ChSaved));
    t2 = tic;

        
for iCh = 1:length(ChSaved)
            
            % Import and crop to block positions
            
            try
                [data_, ~, ~] = load_open_ephys_data(fname_in{iCh});
                data_(1:BlkLength*(iBlk-1))=[]; % Strip off previously processed blocks
%                 
%                 data_(data_>75)=75;
%                 data_(data_<-75)=-75;
                if iBlk ~= noBlks 
                    data{iCh} = data_(1:BlkLength);
                else
                    data{iCh} = data_(1:end);
                end
               
                data_=[];
                flag_{iCh} = 1;   

            catch                 
                fprintf('WARNING Conversion error for channel %s! Padding with zeros.\n',num2str(iCh))
                %                 % default to zero in data(:,iCh)
                flag_{iCh} = 0;   
            end
        end
        data = cell2mat(data);
        info.flag(:,iBlk)= cell2mat(flag_);
        fprintf('>>> Import took %ds seconds, now writing to .bin ...\n',round(toc(t2))) 

        % Write this block to .bin
        t3 = tic;
        fid   = fopen(fname_out{iBlk},'w');
        fwrite(fid,double(data)','int16');
        fclose(fid);
        disp(['>>> Writing block ' num2str(iBlk) ' to .bin file (' fname_out{iBlk} ') took ' sprintf('%ds.',round(toc(t3)))])
        fprintf('>>> Block %d done. Took %ds.\n', iBlk, round(toc(t1))) 
        data = [];
    end  
    
    info.input_directory    = input_directory;
    info.output_directory   = output_directory;
    info.fname_in           = fname_in;
    info.fname_out          = fname_out;
    info.fLength            = fLength;
    info.noBlks             = noBlks;
    info.BlkLength          = BlkLength;
    info.BlkLengthMin       = BlkLengthMin;
    info.BlksizeGB          = BlksizeGB;
	info.ChSaved            = ChSaved;

%finish by concatenating the files.

nBinFiles = numel(dir([PathName '\*.bin']));
binFiles = dir([PathName '\*.bin']);
fileOut = fopen([PathName '\dataALL.bin'],'w');

for i = 1:nBinFiles
    disp(['Writing file Raw_Part_' num2str(i) '.bin'])
    fileIn= fopen([PathName, '/' 'Raw_Part_' num2str(i) '.bin']);
    temp = fread(fileIn);
    fwrite(fileOut,temp);
    fclose(fileIn);
end

fprintf('\n Done')
fclose(fileOut);
