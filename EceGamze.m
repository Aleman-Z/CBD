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

% Change the determined mean interval according to desired one.
desired_mean_interval = 1800;  
%How many secs to read a cell
persec_array_length=21600/length(ratio_sw); %Determined the sampling period
%persec_array_length=round(10800/length(ratio_sw)); %Use this one for 3 hours data
mean_length = round(desired_mean_interval / persec_array_length);
divide_by_sw = round((length(ratio_sw)/mean_length));
difference = divide_by_sw*mean_length - length(ratio_sw);
if difference < 0
    data_manipulated_sw = ratio_sw(1:length(ratio_sw) + difference);
    data_manipulated_theta =  ratio_theta(1:length(ratio_sw) + difference);
else
    data_manipulated_sw = [ratio_sw, zeros(1,difference)];
    data_manipulated_theta =  [ratio_theta, zeros(1,difference)];
end;

data_manipulated_sw = reshape(data_manipulated_sw, divide_by_sw, mean_length);
data_manipulated_theta = reshape(data_manipulated_theta, divide_by_sw, mean_length);
mean_values_sw = mean(data_manipulated_sw, 2);
mean_values_theta = mean(data_manipulated_theta, 2);

time_values = desired_mean_interval*(1:1:divide_by_sw);

x = time_values;
figure
plot (x, mean_values_sw, 'color',[0 0.4470 0.7410]);
hold on
plot(x, mean_values_theta, 'color',[0.8500 0.3250 0.0980]);

hold off
legend('theta','slow wave');
xlabel ('time in secs');
ylabel ('frequency');

% Change the figure name according to the desired_mean_interval
savefig('Rat5_30mins.fig')
    
    