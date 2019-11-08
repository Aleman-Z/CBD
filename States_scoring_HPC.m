%clc
%clear

SampleRate = 30000;  % EEG sampling rate
DS_fac = 100;        % downsampling factor, adjust accordingly

% importing eeg signals
% EEG = abfload('C:\Users\mmlab\Mojtaba\Data\20180926\LFP\2018_09_26_0003.abf','start',1,'stop','e');
LFP_HPC = ans;

% filtering HPC 
HPC_filt_ds = downsample(LFP_HPC, DS_fac);
SampleRate_ds = SampleRate/DS_fac;
[z,p,k] = cheby1(4,0.3,2*[.2 100]/SampleRate,'bandpass');
[sos,g] = zp2sos(z,p,k);
HPC_filt = filtfilt(sos,g,double(HPC_filt_ds));
% HPC_filt_ds = downsample(HPC_filt, DS_fac);
% SampleRate_ds = SampleRate/DS_fac;

% Calculating P_theta/P_T
% [Pratio, start_t, finish_t] = Power_ratio_MT_MovingWin(signal, flo, fhi, win_size, ol_size, nw, SampleRate)
[ratio_theta, ~, start_t, end_t] = Power_ratio_MT_MovingWin(HPC_filt_ds, 3.5, 4.5, 6, 3, 4, SampleRate_ds);
% calcualting P_SW/P_T
[ratio_sw, ~, start_t, end_t] = Power_ratio_MT_MovingWin(HPC_filt_ds, .3, 1.5, 6, 3, 3, SampleRate_ds);

% Calculating means for index_nrem and index_rem
M_ratio_theta = mean(ratio_theta);
M_ratio_sw = mean(ratio_sw);

% scoring rem, nrem and transitory epochs
index_rem = ratio_theta > M_ratio_theta;   % adjust thresholds if necessary
index_nrem = ratio_sw > M_ratio_sw;

index_trans = 1 - index_rem - index_nrem;
index_trans(index_trans==-1) = 1;

index_rem(index_trans==1) = 0;
index_nrem(index_trans==1) = 0;

% Plot figure for ratio SW and theta
plot(ratio_theta,'color',[0 0.4470 0.7410]);
hold on;
plot(ratio_sw,'color',[0.8500 0.3250 0.0980]);
hold off

legend('theta','slow wave');
