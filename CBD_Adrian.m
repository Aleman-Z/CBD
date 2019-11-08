
% selpath = uigetdir('C:\','Select folder with rat ephys data')
% cd(selpath)

% % Self-made downsampling attempt.
% 
% fs=30000; %Ephys sampling frequency
% dname2=uigetdir([],strcat('Select folder where downsampled data should be saved'));
% 
% %Iterations.
% for k=1:length(rats) %Number of rodents.
% 
% %Main folder    
% cd(selpath)
% 
% %Go to rat's folder
% cd(strcat('rat',num2str(rats(k))));
% g=getfolder;
% cd(g{1})
% 
% %Get channels.
% % chan=eval(strcat('channels.Rat',num2str(rats(k))));
% vr=getfield(channels,strcat('Rat',num2str(rats(k))));%Electrode channels. 
% 
% %Find and load ephys channels.
% cfold=dir;
% cfold={cfold.name};
% cfold=cfold(cellfun(@(x) contains(x,'CH'),cfold));
%    
% %HPC
% cf1=cfold(cellfun(@(x) contains(x,num2str(vr(1))),cfold));
% %PFC
% cf2=cfold(cellfun(@(x) contains(x,num2str(vr(2))),cfold));
% 
% if size(cf1,1)~=1 || size(cf2,1)~=1 
%     error('Ambiguous channel')
%     xo
% end
% 
% %Downsampling to 1Khz
% Wn=[500/(fs/2) ]; % Cutoff=500 Hz
% [b,a] = butter(3,Wn); %Filter coefficients for LPF
% 
% %Hippocampus
% [HPC, ~, ~] = load_open_ephys_data_faster(cf1{1});    
% HPC=filtfilt(b,a,HPC);
% HPC=decimator(HPC,fs/1000);
% 
% %PFC
% [PFC, ~, ~] = load_open_ephys_data_faster(cf2{1});
% PFC=filtfilt(b,a,PFC);
% PFC=decimator(PFC,fs/fs_new);
% 
% 
% 
% cd(dname2)
% % if ~exist(num2str(Rat))
% if ~isfolder(num2str(rats(k)))
%     mkdir(num2str(rats(k)))
% end
% cd(num2str(rats(k)))
% 
% 
% if ~exist(str2{num}, 'dir')
%    mkdir(str2{num})
% end
% 
% 
% 
% cd(str2{num})
% if iter_no_saving~=1
%  save('PFC.mat','PFC')
%  save('HPC.mat','HPC')
% %  save('sos.mat','sos')
%  ftext = fopen( str1{num}, 'w' );  
%  fclose(ftext);
% end
% clear PFC HPC %sos
% 
% end
%%

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
vr=getfield(channels,strcat('Rat',num2str(rats(k))));%Electrode channels. 


% %Find and load .mat files.
cfold=dir;
cfold={cfold.name};
cfold=cfold(cellfun(@(x) contains(x,'.mat'),cfold));

%HPC
cf1=cfold(cellfun(@(x) contains(x,num2str(vr(1))),cfold));
'Loading HPC'
load(cf1{1})

% importing eeg signals
LFP_HPC = ans;

% filtering and downsampling HPC 
%Adrian's note: Wrong order of decimation.
%Should be first LPF and then decimation.
fs=30000;
fs_new=300;

Wn=[fs_new/fs ]; % Cutoff=300 Hz
[b,a] = butter(3,Wn); %Filter coefficients for LPF

HPC_filt=filtfilt(b,a,LFP_HPC);
HPC_filt=decimator(HPC_filt,fs/fs_new);

% HPC_filt_ds = downsample(LFP_HPC, DS_fac);
% SampleRate_ds = SampleRate/DS_fac;
% [z,p,k] = cheby1(4,0.3,2*[.2 100]/SampleRate,'bandpass');
% [sos,g] = zp2sos(z,p,k);
% HPC_filt = filtfilt(sos,g,double(HPC_filt_ds));



% Calculating P_theta/P_T
% [Pratio, start_t, finish_t] = Power_ratio_MT_MovingWin(signal, flo, fhi, win_size, ol_size, nw, SampleRate)
[ratio_theta, ~, start_t, end_t] = Power_ratio_MT_MovingWin(HPC_filt_ds, 3.5, 4.5, 6, 3, 4, SampleRate_ds);
% calcualting P_SW/P_T
[ratio_sw, ~, start_t, end_t] = Power_ratio_MT_MovingWin(HPC_filt_ds, .3, 1.5, 6, 3, 3, SampleRate_ds);

% Calculating means for index_nrem and index_rem
M_ratio_theta = mean(ratio_theta);
M_ratio_sw = mean(ratio_sw);




xo
%PFC
cf2=cfold(cellfun(@(x) contains(x,num2str(vr(2))),cfold));

end
