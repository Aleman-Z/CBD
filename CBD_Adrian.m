%Add folder to path and load experiment data. 

% selpath = uigetdir('C:\','Add Adritools folder to path');
% addpath(genpath(selpath))
find_duration=0;
otherstuff=0;
downsample_data=1;

selpath = uigetdir('C:\','Add CBD github folder to path');
addpath(genpath(selpath))

cd(selpath)
load('CBD_2.mat')
clear selpath

%Go to folder with CBD Matlab files.
selpath = uigetdir('C:\','Select CBD_acutes folder');
cd(selpath)


%Starting times for Rat2-Rat5.
t1(1) = duration([12 0 0]);
t1(2) = duration([11 53 0]);
t1(3) = duration([12 0 0]);
t1(4) = duration([11 54 0]);

t2=t1+seconds(durations_trials)'; %Values calculated from a previous iteration.

if find_duration==1
durations_trials=zeros(length(rats),1);
end

%Iterations.
for r=1:length(rats) %Number of rodents.

%Main folder    
cd(selpath)

%Go to rat's folder
cd(strcat('Rat',num2str(rats(r))));

%Get channels.
% chan=eval(strcat('channels.Rat',num2str(rats(k))));
chan=getfield(channels,strcat('Rat',num2str(rats(r))));%Electrode channels. 


% %Find and load .mat files.
cfold=dir;
cfold={cfold.name};
cfold=cfold(cellfun(@(x) ~isempty(strfind(x,'.mat')),cfold));

%Parameters
fs=30000; %Sampling frequency of acquisition.
fs_new=1000; %New sampling frequency. Max freq= fs_new/2.

%HPC
cf1=cfold(cellfun(@(x) ~isempty(strfind(x,num2str(chan(1)))),cfold));
'Loading HPC'
load(cf1{1})

clear ans HPC_filt HPC_filt_ds

% % importing eeg signals
% LFP_HPC = ans;

%For duration uncomment this:
if find_duration==1
    durations_trials(r,:)=length(LFP_HPC)/fs;
    continue 
end

% filtering and downsampling HPC 
%Adrian's note: Wrong order of decimation.
%Should be first filtering and then decimation.

%Butter filter gives smoother signal than Chebyshev.
Wn=[fs_new/fs ]; % Cutoff=500 Hz
[b,a] = butter(3,Wn); %Filter coefficients for LPF

HPC_filt=filtfilt(b,a,LFP_HPC);
% HPC_filt=decimator(HPC_filt,fs/fs_new);
% HPC_filt=HPC_filt';
clear LFP_HPC
HPC_downsampled=downsample(HPC_filt,fs/fs_new);
clear HPC_filt

if downsample_data==1
    save('HPC_downsampled2.mat','HPC_downsampled')
    
%PFC
cf1=cfold(cellfun(@(x) ~isempty(strfind(x,num2str(chan(2)))),cfold));
'Loading PFC'
load(cf1{1})
clear ans HPC_filt HPC_filt_ds

% importing eeg signals
LFP_PFC = LFP_HPC;
clear LFP_HPC

PFC_filt=filtfilt(b,a,LFP_PFC);
clear LFP_PFC
PFC_downsampled=downsample(PFC_filt,fs/fs_new);
clear PFC_filt

save('PFC_downsampled2.mat','PFC_downsampled')

continue

end

% HPC_filt_ds = downsample(LFP_HPC, DS_fac);
% SampleRate_ds = SampleRate/DS_fac;
% [z,p,k] = cheby1(4,0.3,2*[.2 100]/SampleRate,'bandpass');
% [sos,g] = zp2sos(z,p,k);
% HPC_filt = filtfilt(sos,g,double(HPC_filt_ds));

nw=1.25; %Minimum taper value.
%Too many tapers smooth the signal. 
% 6 sec window, 3 sec overlap (50%)
% Calculating P_theta/P_T
% [Pratio, start_t, finish_t] = Power_ratio_MT_MovingWin(signal, flo, fhi, win_size, ol_size, nw, SampleRate)
[ratio_theta, ~, ~ , ~] = Power_ratio_MT_MovingWin(HPC_filt, 3.5, 4.5, 6, 3, nw, fs_new);
% calcualting P_SW/P_T
[ratio_sw, ~, ~,~] = Power_ratio_MT_MovingWin(HPC_filt, .3, 1.5, 6, 3, nw, fs_new);

M=[ratio_theta; ratio_sw];
save('M.mat','M')
clear M ratio_sw ratio_theta 

if otherstuff==1
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
    %Here 21600 is the equivalent of 6 hours in seconds. 
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

    % time_values = desired_mean_interval*(1:1:divide_by_sw);
    % 
    % x = time_values;
    % figure
    % plot (x, mean_values_sw, 'color',[0 0.4470 0.7410]);
    % hold on
    % plot(x, mean_values_theta, 'color',[0.8500 0.3250 0.0980]);
    % 
    % hold off
    % legend('theta','slow wave');
    % xlabel ('time in secs');
    % ylabel ('frequency');

    % Change the figure name according to the desired_mean_interval
    % savefig('Rat5_30mins.fig')
    % %HPC Normalized values
    % N_ratio_theta = (ratio_theta-mean(ratio_theta))/std(ratio_theta);
    % N_ratio_sw =(ratio_sw-mean(ratio_sw))/std(ratio_sw);
    % HPC=[N_ratio_theta; N_ratio_sw];


    %PFC
    % cf2=cfold(cellfun(@(x) contains(x,num2str(vr(2))),cfold));
end

end
