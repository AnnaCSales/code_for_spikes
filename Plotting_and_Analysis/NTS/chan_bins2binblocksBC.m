function info = OE2binblocks(varargin)
% Aleks Domanski UoB 2016 hacked away at a bit by Anna Sales 2017
% Parallel import and chunking of Open Ephys data into equal length .bin
% files. NB the last .bin in the sequence will appear smaller since it
% contains the shorter block of data at the end of the recording (i.e .the remainder) 

% Corrected 03032018 to account for >9 bin files - now concatenates
% reliably in the right order. Anna Sales.

addpath(genpath('/panfs/panasas01/phph/phacs/MATLAB/analysis-tools-master')); % path to OEP utilities folder
addpath(genpath('/panfs/panasas01/phph/phacs/MATLAB/npy-matlab-master')); % path to npy-matlab scripts
addpath(genpath('/panfs/panasas01/phph/phacs/MATLAB')); % path to my own matlab scripts
 
folder_name='';  %Update this before running
[PathName] = pwd;% ['/panfs/panasas01/phph/phacs/Data/20180221/' folder_name];  
input_directory=PathName;
output_directory = input_directory;

% Specify the bin file length in minutes here:
BlkLengthMin  = 10; 


info = get_session_info(input_directory);


%% Brute force channel import 
    %N.B. 10 minutes of 32 channel recording @ 16bit/30kHz  ~ 1.15GB
      
ChSaved = 1:32;
processor_index = 0;

% grab a channel to get total recording length
fname_in{1} = [input_directory filesep '100_CH1.continuous'];  %CAN AMEND THIS FOR RECORDINGS WITH _2 etc
evalc('[~,~, fileInfo] = load_open_ephys_data(fname_in{1});'); % Supress commandline output w/ evalc

fLength       = fileInfo.header.bufferSize * length (fileInfo.ts);   % File length in samples
BlkLength     = BlkLengthMin * 60 * fileInfo.header.sampleRate;      % Block length in samples
BlksizeGB     = (BlkLength * 16*  numel(ChSaved))/(1e9*8);           % Block size in GB
noBlks        = ceil(fLength/BlkLength);                             % Number of data blocks to process

fprintf('>>> Recording is %d seconds long (%.2f hours).\n',round(fLength/fileInfo.header.sampleRate),round(fLength/fileInfo.header.sampleRate)/3600 )
fprintf('>>> Writing data as %d blocks of %d minutes, each block is %.3f GB.\n',noBlks,BlkLengthMin,BlksizeGB)

% Cache filenames
for iCh =  1:length(ChSaved)
    fname_in{iCh} = [input_directory filesep '100_CH' int2str(ChSaved(iCh)) '.continuous'];
end
for iBlk = 1:noBlks % deal with full blocks first
    fname_out{iBlk} = [get_full_path(output_directory) filesep 'Raw_Part_' num2str(iBlk) '.bin'];
end

% Loop across time in blocks
for iBlk = 1:noBlks
   
    fprintf('>>> Starting Block %d of %d.\n', iBlk, ceil(fLength/BlkLength))

    % Preallocate block data (special case for last block as will be time remainder)
    if iBlk ~= noBlks
        for iCh = 1 : length(ChSaved)
            data{iCh} = zeros(BlkLength,1);

        end
        %             data = zeros(BlkLength,numel(ChSaved));
        %             data = mat2cell(data,size(data,1),ones(1,size(data,2)));
    else
        for iCh = 1 : length(ChSaved)
            data{iCh} =  zeros(fLength-BlkLength*(noBlks-1),1);
        end
        %            data = zeros(fLength-BlkLength*(noBlks-1),numel(ChSaved));
        %             data = mat2cell(data,size(data,1),ones(1,size(data,2)));
    end
    flag_ = cell(1,length(ChSaved));
    t2 = tic;


for iCh = 1:length(ChSaved)
            
            % Import and crop to block positions
            
            try
                [data_, ~, ~] = load_open_ephys_data(fname_in{iCh});
                data_(1:BlkLength*(iBlk-1))=[]; % Strip off previously processed blocks
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

nBinFiles = numel(dir([PathName '/*.bin']));
binFiles = dir([PathName '/*.bin']);
fileOut = fopen([PathName, '/100_raw.bin'],'w');

for i = 1:nBinFiles
    disp(['Writing file Raw_Part_' num2str(i) '.bin'])
    fileIn= fopen([PathName, '/' 'Raw_Part_' num2str(i) '.bin']);
    temp = fread(fileIn);
    fwrite(fileOut,temp);
    fclose(fileIn);
end

fclose(fileOut);


 delete(gcp)