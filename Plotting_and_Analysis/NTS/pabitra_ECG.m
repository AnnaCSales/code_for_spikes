dataCAR=load_open_ephys_binary('structure.oebin', 'continuous',2,'mmap'); % 1 is FGPA, 2 is CAR
fs=dataCAR.Header.sample_rate;
bitVolts_probe=dataCAR.Header.channels(1).bit_volts;
tbase=double(dataCAR.Timestamps) / fs;

dataRAW=load_open_ephys_binary('structure.oebin', 'continuous',1,'mmap'); % 1 is FGPA, 2 is CAR
fs=dataRAW.Header.sample_rate;
bitVolts_ADC=dataRAW.Header.channels(34).bit_volts;
ADC2=dataRAW.Data.Data.mapped(34,1:end);
raw_signal=double(ADC2)*bitVolts_ADC;
 
LPlim=40;   %filter below 40Hz
[Db,Da]=butter(2,LPlim/(0.5*fs) );
LP_ECG=filtfilt(Db,Da,raw_signal);  %filtered delta
resp=medfilt1(raw_signal);

%%
 
plotwin=[300, 310];
plotind=find(tbase>plotwin(1) & tbase< plotwin(2));

figure('Color', 'w', 'Units', 'Normalized', 'Position', [0.1 0.1 0.5 0.75])
hold on
a=subplot(2,1,1);

hold on
for c=[1,5,13,25]
    plot_level=(c-1)*100;
    this_data=double(dataCAR.Data.Data.mapped(c, plotind)) * bitVolts_probe;
    plot(tbase(plotind), plot_level+this_data, 'b');
    
end

xlabel('Time (s)')
ylabel('Probe data')
title('Probe data', 'Fontweight', 'normal')
% ylim([-1000 16000])

b=subplot(2,1,2);
plot(tbase(plotind), resp(plotind))
box off
xlabel('Time (s)')
ylabel('\muV')
title('ECG', 'Fontweight', 'normal')
ylim([-2.5, 2])
a.Position=[0.1 0.38 0.8 0.55];
b.Position=[0.1 0.06 0.8 0.2];
%%
plotwin=[300, 310];
plotind=find(tbase>plotwin(1) & tbase< plotwin(2));

figure('Color', 'w', 'Units', 'Normalized', 'Position', [0.1 0.1 0.5 0.75])
hold on
a=subplot(4,1,1);

hold on
plt=1;
for c=[5,13,25]
    subplot(4,1,plt)
%     plot_level=(c-1)*500;
    this_data=double(dataCAR.Data.Data.mapped(c, plotind)) * bitVolts_probe;
    plot(tbase(plotind),this_data, 'b');
    title(['CH ' num2str(c)], 'FontWeight', 'normal')
    box off
    ylabel('\muV')
    plt=plt+1;
end


ylabel('\muV')
% ylim([-1000 16000])

subplot(4,1,4);
plot(tbase(plotind), resp(plotind))
box off
xlabel('Time (s)')
ylabel('\muV')
title('ECG /resp', 'Fontweight', 'normal')
ylim([-2.5, 2])

%%
ADC2=dataRAW.Data.Data.mapped(34,1:end);
raw_signal=double(ADC2)*bitVolts_ADC;
 
LPlim=40;   %filter below 40Hz
[Db,Da]=butter(2,LPlim/(0.5*fs) );

LP_ECG=filtfilt(Db,Da,raw_signal);  %filtered delta

subplot(2,1,2)
plot(tbase(plotind), LP_ECG(plotind))

resp=medfilt1(raw_signal);

figure('Color', 'w')
subplot(2,1,1)
plot(tbase(plotind), LP_ECG(plotind))
subplot(2,1,2)
plot(tbase(plotind), resp(plotind))

