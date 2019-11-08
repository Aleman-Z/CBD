clc    %clean the command window
clear  %remove all variables

% hippocampal pyr layer. Theta
SampleRate = 30000;  % EEG sampling rate
DS_fac = 300;        % downsampling factor, adjust accordingly
[LFP_data, ~, ~] = load_open_ephys_data('100_CH4_0.continuous');  % include the channel according to the animal and the cortical area
LFP_data_r = resample(LFP_data, 1, DS_fac);

% filtering HPC 
[z,p,k] = cheby1(4,0.3,2*[.2 100]/SampleRate,'bandpass');
[sos,g] = zp2sos(z,p,k);
HPC_filt = filtfilt(sos,g,double(LFP_data));
HPC_filt_ds = downsample(HPC_filt, DS_fac);
SampleRate_ds = SampleRate/DS_fac;

% Calculating P_theta/P_T
 [ratio_theta, ~, start_t, end_t] = Power_ratio_MT_MovingWin(HPC_filt_ds, 3.5, 4.5, 6, 3, 4, SampleRate_ds);
 
% calcualting P_SW/P_T
[ratio_sw, ~, start_t, end_t] = Power_ratio_MT_MovingWin(HPC_filt_ds, .3, 1.5, 6, 3, 3, SampleRate_ds);

% scoring rem, nrem and transitory epochs
index_rem = ratio_theta > 0.4;   % adjust thresholds if necessary
index_nrem = ratio_sw > 0.4;

%normalized theta ratio
mean_theta = mean(ratio_theta); %calculating the mean of theta
theta_nor = ratio_theta - mean_theta; %theta ratio normalized

%normalized delta ratio
mean_sw = mean(ratio_sw); %calculating the mean of theta
sw_nor = ratio_sw - mean_sw; %sw ratio normalized
%Plotting 
plot(theta_nor, 'k-')
hold on
plot(sw_nor, 'c-')
savefig('theta_sw')

