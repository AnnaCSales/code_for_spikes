% %% Takes 40 chan binary and reformats into32, for when you accidentally record the ADC
addpath('D:\Code\MATLAB\Kilosort_and_binmaking')

dataRAW=load_open_ephys_binary('structure.oebin', 'continuous',1,'mmap'); % 1 is FGPA, 2 is CAR
fs=dataRAW.Header.sample_rate;
bitVolts_ADC=dataRAW.Header.channels(34).bit_volts;
bitVolts_probe=dataRAW.Header.channels(1).bit_volts;
tbase=double(dataRAW.Timestamps) / fs;
%% Make separate files for each real probe chan

for c=1:32
    fprintf('\n ..Working on channel %d' ,c)
    fn=['channel_' num2str(c)];
    this_data=double(dataRAW.Data.Data.mapped(c, :));
    fileID = fopen([pwd '\raw_data_1to32\' fn],'w');
    fwrite(fileID,this_data,'int16');
    fclose(fileID);
end
%% Sanity check - the two columns should be identical
% 
fid=fopen([pwd '\raw_data_1to32\' 'channel_5'],'r');
ch1=fread(fid, 'int16');
[dataRAW.Data.Data.mapped(1,1:10)', ch1(1:10)]

fid=fopen([pwd '\raw_data_1to32\' 'channel_15'],'r');
ch1=fread(fid, 'int16');
[dataRAW.Data.Data.mapped(15,1:10)', ch1(1:10)]
%% Now make bin

folder_name='';  %Update this before running
PathName =[ pwd '\raw_data_1to32\'];% 
input_directory=PathName;
output_directory = input_directory;

% Specify the bin file length in minutes here:
BlkLengthMin  = 10; 
ChSaved = 1:32;
fileID = fopen([PathName 'channel_1']);
testfile = fread(fileID, 'int16');
fclose(fileID);
    
fLength       = numel(testfile);   % File length in samples
BlkLength     = BlkLengthMin * 60 * fs;      % Block length in samples
BlksizeGB     = (BlkLength * 16*  numel(ChSaved))/(1e9*8);           % Block size in GB
noBlks        = ceil(fLength/BlkLength);                             % Number of data blocks to process
  
fprintf('>>> Recording is %d seconds long (%.2f hours).\n',round(fLength/fs),round(fLength/fs/3600) )
fprintf('>>> Writing data as %d blocks of %d minutes, each block is %.3f GB.\n',noBlks,BlkLengthMin,BlksizeGB)

% Cache filenames
for iCh =  1:length(ChSaved)
    fname_in{iCh} = [pwd '\raw_data_1to32\channel_' num2str(iCh)];
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
            fname_in{iCh}
            fileID_ = fopen(fname_in{iCh});
            data_ = fread(fileID, 'int16');
            fclose(fileID_);
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
    fwrite(fid, data','int16');
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

%% Sanity check, again
dataRAW.Data.Data.mapped(1,1:10)'

m=memmapfile( [pwd '\raw_data_1to32\100_raw.bin'],'Format','int16') ;
ch1_new=m.Data(1:32:3300)';

[dataRAW.Data.Data.mapped(1,1:10)', ch1_new(1:10)']
%%
applyCARtoDatAndFilter('100_raw.bin',32, pwd);