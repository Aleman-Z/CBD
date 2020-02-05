clear
clc
% %Open Ephys
% addpath(uigetdir('C:\','Select Open Ephys functions folder'));

%Fieldtrip
addpath('C:\Users\students\Documents\Tugdual\GitHub\fieldtrip');

%CorticoHippocampal
addpath('C:\Users\students\Documents\Tugdual\GitHub\CorticoHippocampal');

%ADRITOOLS
addpath('C:\Users\students\Documents\Tugdual\GitHub\ADRITOOLS');


% patID=uigetdir('C:\','Select CBD folder');
patID = 'C:\Users\students\Documents\Tugdual\GitHub\CBD'
cd(patID)
% cd('CBD')
load('CBD_5.mat')
%% SELECT RAT
opts.Resize = 'on';
opts.WindowStyle = 'modal';
opts.Interpreter = 'tex';
prompt=strcat('\bf Select a rat#. Options:','{ }',num2str(rats));
answer = inputdlg(prompt,'Input',[2 30],{''},opts);
Rat=str2double(answer{1});

%% SELECT BRAIN AREA
BA = inputdlg({'Brain area from which you would like the power spectrum'},...
              'Type your selection', [1 60]);
          
%% LOAD DATA

cfold=dir(strcat('C:\Users\students\Documents\Tugdual\Downsampled_data\',answer{1},'\CBD'));
cfold={cfold.name};
cfold=cfold(cellfun(@(x) contains(x,BA{1}),cfold));


Data = load(strcat('C:\Users\students\Documents\Tugdual\Downsampled_data\',answer{1},'\CBD\',cfold{1}));
Data = getfield(Data, BA{1});

%% DIVIDE INTO 1s EPOCHS

fs = 600; %Hz
e_t = 1; %s
e_samples = e_t*(fs); 
 
ch = length(Data);

nc = floor(ch/e_samples); %Number of epochs
NC=[];
for kk=1:nc
    NC(:,kk)= Data(1+e_samples*(kk-1):e_samples*kk);
end

max_nc = [];
for k=1:size(NC,2)
%      plot(NC(:,k))
%      hold on
    max_nc(k) = max(NC(:,k));
end

outliers = isoutlier(max_nc, 'quartiles');

prefiltered_NC = NC(:, ~outliers);


%% POWER SPECTRUMS

[pxx,f] = pmtm(prefiltered_NC,4,[],fs);
px=mean(pxx,2);
s = semilogy(f,(px).*f,'Color',[0 0 0],'LineWidth',2);
grid on
title(strcat('Power spectrum of the', {' '}, BA{1},' of rat',{' '}, answer{1}))
xlabel('Frequency (Hz)')
ylabel('Power')

saveas(s,strcat('C:\Users\students\Documents\Tugdual\PowerSpectrums\',BA{1},answer{1},'.png'))