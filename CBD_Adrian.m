
%Add folder to path and load experiment data. 

% selpath = uigetdir('C:\','Add Adritools folder to path');
% addpath(genpath(selpath))

selpath = uigetdir('C:\','Add CBD github folder to path');
addpath(genpath(selpath))

cd(selpath)
load('CBD_1.mat')
clear selpath

%Go to folder with CBD Matlab files.
selpath = uigetdir('C:\','Select CBD_acutes folder');
cd(selpath)

%Iterations.
for k=1:length(rats) %Number of rodents.

%Main folder    
cd(selpath)

%Go to rat's folder
cd(strcat('Rat',num2str(rats(k))));

%Get channels.
% chan=eval(strcat('channels.Rat',num2str(rats(k))));
chan=getfield(channels,strcat('Rat',num2str(rats(k))));%Electrode channels. 


% %Find and load .mat files.
cfold=dir;
cfold={cfold.name};
cfold=cfold(cellfun(@(x) ~isempty(strfind(x,'.mat')),cfold));

%Parameters
fs=30000; %Sampling frequency of acquisition.
fs_new=300; %New sampling frequency. Max freq= fs_new/2.

%HPC
cf1=cfold(cellfun(@(x) ~isempty(strfind(x,num2str(chan(1)))),cfold));
'Loading HPC'
load(cf1{1})

% importing eeg signals
LFP_HPC = ans;

% filtering and downsampling HPC 
%Adrian's note: Wrong order of decimation.
%Should be first filtering and then decimation.

%Butter filter gives smoother signal than Chebyshev.
Wn=[fs_new/fs ]; % Cutoff=300 Hz
[b,a] = butter(3,Wn); %Filter coefficients for LPF

HPC_filt=filtfilt(b,a,LFP_HPC);
% HPC_filt=decimator(HPC_filt,fs/fs_new);
% HPC_filt=HPC_filt';
HPC_filt=downsample(HPC_filt,fs/fs_new);

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

% % Calculating means for index_nrem and index_rem
% M_ratio_theta = mean(ratio_theta);
% M_ratio_sw = mean(ratio_sw);

%HPC Normalized values
N_ratio_theta = (ratio_theta-mean(ratio_theta))/std(ratio_theta);
N_ratio_sw =(ratio_sw-mean(ratio_sw))/std(ratio_sw);
HPC=[N_ratio_theta; N_ratio_sw];


%PFC
cf2=cfold(cellfun(@(x) contains(x,num2str(vr(2))),cfold));

end
