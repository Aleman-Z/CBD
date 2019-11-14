% S1=Matrix{1};
% S2=Matrix{2};
% %%
% plot(linspace(t1(1),t2(1),length(S1(1,:))),S1(1,:))
% hold on
% plot(linspace(t1(2),t2(2),length(S2(1,:))),S2(2,:))
% 
% %%
% imagesc(S1(:,:))
%%
selpath = uigetdir('C:\','Add CBD github folder to path');
addpath(genpath(selpath))

cd(selpath)
load('CBD_2.mat')
clear selpath

selpath = uigetdir('C:\','Select folder with Downsampled_data');

plot_allanimals=1;
%Starting times for Rat2-Rat5.
t1(1) = duration([12 0 0]);
t1(2) = duration([11 53 0]);
t1(3) = duration([12 0 0]);
t1(4) = duration([11 54 0]);

t2=t1+seconds(durations_trials)'; %Values calculated from a previous iteration.
% xo
freq_band.RawSignal=[0 150-1];
freq_band.Delta=[0.1 4];
freq_band.Theta=[4 8];
freq_band.SpindlesRange=[10 20];
freq_band.HighGamma=[80 150-1]; %250 
FN=fieldnames(freq_band);

fs=300; %Sampling freq after downsampling.
%% 
close all
for lab=1:length(label1)
    for r=1:length(rats)
        for k=1:length(fieldnames(freq_band))

        % cd('C:\Users\addri\Documents\Donders\Projects\CBD\Downsampled_data')
        cd(selpath)

        % %Find and access proper folder.
        cfold=dir;
        cfold={cfold.name};
        cf1=cfold(cellfun(@(x) ~isempty(strfind(x,num2str(rats(r)))),cfold));
        cd(cf1{1})

        % %Find and load .mat files.
        cfold=dir;
        cfold={cfold.name};
        cfold=cfold(cellfun(@(x) ~isempty(strfind(x,'.mat')),cfold));
        cfold=cfold(cellfun(@(x) ~isempty(strfind(x,label1{lab})),cfold));
        signal_ds=load(cfold{1});
        signal_ds=getfield(signal_ds,erase(cfold{1},'.mat'));

        %Bandpassing.
        if ~strcmp(FN{k},'RawSignal')
            GF=getfield(freq_band,FN{k});
            Wn=[GF(1)/(fs/2) GF(2)/(fs/2) ]; % Sampling freq=300 Hz.
            [b,a] = butter(3,Wn); %Filter coefficients for LPF
            signal_ds=filtfilt(b,a,signal_ds);    
        end
        
        signal_normal=(signal_ds-mean(signal_ds))/std(signal_ds);

            if plot_allanimals==1;
                plot(linspace(t1(r),t2(r),length(signal_normal)),signal_normal+100*k)
                hold on
                xlabel('Time')
            else
                %% Epoching data
                T1=1;
                T2=T1+2*(300);
                con=1;
                while T2 <= length(signal_normal)
                %T2=T1+2*(300);    
                  new_seg = signal_normal(T1:T2);
                  NC(con,:)=new_seg;

                con=con+1;
                T2=T2+2*(300);
                T1=T1+2*(300);
                end
            [pxx,f]= periodogram(NC.',hann(size(NC.',1)),size(NC.',1),300);    
            %%
            px=mean(pxx,2);
            s=semilogy(f,(px),'LineWidth',2);
            xlabel('Frequency (Hz)')
            ylabel('Power')
            printing(strcat(label1{lab},'_periodogram'))

            end
        % ylabel('Rat')
        end
        
        if plot_allanimals==1;
            title(strcat('Normalized',{' '},label1{lab},{' '}, 'recording',{' '},FN{k}))
            yticks([100 200 300 400 500])
            yticklabels({FN{1},FN{2},FN{3},FN{4},FN{5}})
            xo
%             saveas(gcf,strcat(label1{lab},'_traces',FN{k},'.fig'))

            close all
        end
    end
end



% end    
    %    xo
%% PFC
% 
% for r=1:length(rats)
% 
% cd(selpath)
% 
% % %Find and access proper folder.
% cfold=dir;
% cfold={cfold.name};
% cf1=cfold(cellfun(@(x) ~isempty(strfind(x,num2str(rats(r)))),cfold));
% cd(cf1{1})
% 
% % %Find and load .mat files.
% cfold=dir;
% cfold={cfold.name};
% cfold=cfold(cellfun(@(x) ~isempty(strfind(x,'.mat')),cfold));
% load(cfold{2}) %PFC
% 
% PFC_normal=(PFC_downsampled-mean(PFC_downsampled))/std(PFC_downsampled);
% 
%     if plot_allanimals==1;
%         plot(linspace(t1(r),t2(r),length(PFC_normal)),PFC_normal+100*r)
%         hold on
%         xlabel('Time')
%         % ylabel('Rat')
%     else
%         %% Epoching data
%         T1=1;
%         T2=T1+2*(300);
%         con=1;
%         while T2 <= length(PFC_normal)
% %         T2=T1+2*(300);    
%           new_seg = PFC_normal(T1:T2);
%           NC(con,:)=new_seg;
% 
%         con=con+1;
%         T2=T2+2*(300);
%         T1=T1+2*(300);
% %         T1=T2;
%         end
%     [pxx,f]= periodogram(NC.',hann(size(NC.',1)),size(NC.',1),300);    
%     %%
%     px=mean(pxx,2);
%     s=semilogy(f,(px),'LineWidth',2);
%     xlabel('Frequency (Hz)')
%     ylabel('Power')
%     printing('PFC_periodogram')
%     
% end
% 
%     end
% %end
% 
%     if plot_allanimals==1;
%         title('Normalized PFC recording')
%         yticks([100 200 300 400])
%         yticklabels({'Rat 2','Rat 3','Rat 4','Rat 5'})
% %         printing('PFC_traces')
%         saveas(gcf,strcat('PFC_traces','.fig'))
%         close all
%     end

% %% Periodogram
% [pxx,w] = periodogram(PFC_normal);
% plot(w,10*log10(pxx))
% 
% %% Epoching data
% T1=1;
% T2=1;
% con=1;
% while T2 <= length(PFC_normal)
% T2=T1+2*(300);    
%   new_seg = PFC_normal(T1:T2);
%   NC(con,:)=new_seg;
% 
% con=con+1;
% T1=T2;
% end

%%

% [pxx,f]= periodogram(NC.',hann(size(NC.',1)),size(NC.',1),300);    
% %%
% px=mean(pxx,2);
% s=semilogy(f,(px),'LineWidth',2);
% xlabel('Frequency (Hz)')
% ylabel('Power')
% % title(strcat('Power in NREM',{' '} ,label1{2*w-1} ,{' '},'signals'))