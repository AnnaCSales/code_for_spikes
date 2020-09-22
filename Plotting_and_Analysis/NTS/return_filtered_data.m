function [ecg_sig, resp_sig] = return_filtered_data(raw_sig, fs)
%Takes in raw ADC data from the NTS prep, returns a signal filtered to
%clearly display respiratory activity / ECG respectively.

% Start with the resp. signal:

zsig=zscore(raw_sig); %find the zscored version of the raw signal
mn_sg=mean(zsig); %find the mean of the signal 
std_sig=std(zsig); %find the standard deviation of the signal 

%set anything above 2*stds back to the mean (remove high ECG readings):
nstd=2;
high_ind=find(abs(zsig)>nstd*std_sig);
zsig(high_ind)=0;

%take the moving mean. This will get rid of the higher frequency ECG
%components:
resps=movmean(abs(zsig),8000);

%now high pass filter to remove slow drifts in the baseline:
HPlim=0.2;   %the frequency limit for the HP filter, in Hz:
[Db,Da]=butter(2,HPlim/(0.5*fs), 'high' );
resp_sig=filtfilt(Db,Da,resps);  %Our final processed signal.

% Now filter for ECG 
 
% LP filter to remove laser artefact / clean up for peak
% detection

HPlim=100;   %the frequency limit for the LP filter, in Hz (100):
[Db,Da]=butter(2,HPlim/(0.5*fs), 'low' );
ecg_sig=filtfilt(Db,Da,raw_sig);  %Our final processed signal.

end

