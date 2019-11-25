
%%
selpath = uigetdir('C:\','Add CBD github folder to path');
addpath(genpath(selpath))

cd(selpath)
load('CBD_4.mat')
clear selpath

selpath = uigetdir('C:\','Select folder with Downsampled_data');

plot_allanimals=1;
%Starting times for Rat2-Rat5.
t1(1) = duration([12 0 0]);
t1(2) = duration([11 53 0]);
t1(3) = duration([12 0 0]);
t1(4) = duration([11 54 0]);

%Starting times for Rat9-Rat16.
t1(5) = duration([16 20 0]);
t1(6) = duration([15 21 0]);
t1(7) = duration([15 45 0]);
t1(8) = duration([17 18 0]);
t1(9) = duration([11 30 0]);
t1(10) = duration([11 30 0]);
t1(11) = duration([11 24 0]);
t1(12) = duration([11 35 0]);

t2=t1+seconds(durations_trials)'; %Values calculated from a previous iteration.
% xo
freq_band.RawSignal=[0 500-1];
freq_band.Delta=[0.11 4];
freq_band.Theta=[4 8];
freq_band.SpindlesRange=[10 20];
freq_band.HighGamma=[80 250-1]; %250 
FN=fieldnames(freq_band);

fs_new=500; %Sampling freq after downsampling.
amp_vec=[7 4 4 4 5];

colorvec=[
       24    157    233
         233    59    24
%     233    224    24
% 0  0  0
255 0 255
    125    24         233
    36    233         24

];
colorvec=colorvec./255;

rat_condition={
'CBD'    
'Vehicle'
'Vehicle'
'CBD'
'Vehicle'
'CBD'
'CBD'
'Vehicle'
'Vehicle'
'Vehicle'
'CBD'
'CBD'
    };
%% 
close all
% stats_vec=nan(length(rats),length(label1),length(fieldnames(freq_band)),2); %Declare vector to save values.
% stats_vec_norm=nan(length(rats),length(label1),length(fieldnames(freq_band)),2); %Declare vector to save values.

%for lab=1:length(label1)
    for r=1:length(rats)
        stats_vec=nan(length(label1),length(fieldnames(freq_band)),2); %Declare vector to save values.
        stats_vec_norm=nan(length(label1),length(fieldnames(freq_band)),2); %Declare vector to save values.
       
        for lab=1:length(label1)
        for k=1:length(fieldnames(freq_band))

        % cd('C:\Users\addri\Documents\Donders\Projects\CBD\Downsampled_data')
        cd(selpath)

        % %Find and access proper folder.
        cfold=dir;
        cfold={cfold.name};
        cfold=cfold(cellfun(@(x) ~isempty(strfind(x,strcat('Rat',num2str(rats(r))))),cfold));
        cf1=cfold(cellfun(@(x) isempty(strfind(x,num2str('bands'))),cfold));
        cd(cf1{1})

        % %Find and load .mat files.
        cfold=dir;
        cfold={cfold.name};
        cfold=cfold(cellfun(@(x) ~isempty(strfind(x,'.mat')),cfold));

%         cfold=cfold(cellfun(@(x) ~isempty(strfind(x,label1{lab})),cfold));
%         signal_ds=load(cfold{1});
%         signal_ds=getfield(signal_ds,erase(cfold{1},'.mat'));
        if  ~isempty(cfold(cellfun(@(x) ~isempty(strfind(x,label1{lab})),cfold))) %When both areas were recorded
            cfold=cfold(cellfun(@(x) ~isempty(strfind(x,label1{lab})),cfold));
        end
            signal_ds=load(cfold{1});
            signal_ds=getfield(signal_ds,erase(cfold{1},'.mat'));
        if isempty(cfold(cellfun(@(x) ~isempty(strfind(x,label1{lab})),cfold)))
           signal_ds=signal_ds.*0; % Since signal was not recorded make it zero.
        end

        %Bandpassing.
        if ~strcmp(FN{k},'RawSignal')
            GF=getfield(freq_band,FN{k});
            Wn=[GF(1)/(fs_new/2) GF(2)/(fs_new/2) ]; % Sampling freq=500 Hz.
            [b,a] = butter(3,Wn); %Filter coefficients for LPF
            signal_ds=filtfilt(b,a,signal_ds);    
        end
stats_vec(lab,k,1)=mean(signal_ds);
stats_vec(lab,k,2)=std(signal_ds);

        %Normalized version
       if ~isempty(cfold(cellfun(@(x) ~isempty(strfind(x,label1{lab})),cfold))) %If signal is not zero.
        signal_normal=(signal_ds-mean(signal_ds))/std(signal_ds);
       else
        signal_normal=signal_ds;
       end        
%        signal_normal=(signal_ds-mean(signal_ds))/std(signal_ds);
        clear signal_ds
        uplim=repmat(mean(signal_normal)+2*std(signal_normal),length(signal_normal),1);
        lowlim=repmat(mean(signal_normal)-2*std(signal_normal),length(signal_normal),1);
        
%         signal_normal=(signal_ds);
stats_vec_norm(lab,k,1)=mean(signal_normal);
stats_vec_norm(lab,k,2)=std(signal_normal);
%xo
            if plot_allanimals==1;
                if lab==1
                        if k==1
                            allscreen()% Enlarges figure size to fit screen. Download it from ADRITOOLS github.
                                if r==8  %HPC raw for rat 12.
                                    
                                        fact=2;    
                                        uplim=mean(signal_normal)+fact*std(signal_normal);
                                        lowlim=mean(signal_normal)-fact*std(signal_normal);
                                        %xo
                                        nsig=envelope(signal_normal);
                                        %nsig=(signal_normal>uplim(1)| signal_normal<lowlim(1));
                                        nsig=(nsig>uplim(1)| nsig<lowlim(1));
                                        nonZeroElements = nsig ~= 0;
                                        %minSeparation = 50000;
                                        %minSeparation = 10000;
                                        minSeparation=60*(fs_new); %10 seconds

                                        nonZeroElements = ~bwareaopen(~nonZeroElements, minSeparation);
                                        [labeledRegions, numRegions] = bwlabel(nonZeroElements);
                                        labeledRegions(labeledRegions~=0)=1;

                                        %Minimum of 3 sec duration.
                                        labeledRegions=min_duration(labeledRegions,3,fs_new);

                                        h=area(linspace(t1(r),t2(r),length(signal_normal)),labeledRegions.*1100);
                                        h.FaceColor=[1 0 0];
                                        h.FaceAlpha=0.2;

                                    
                                end
                        end
                        plot(linspace(t1(r),t2(r),length(signal_normal)),signal_normal*amp_vec(k)+100*k,'Color',colorvec(k,:))
                        hold on
                        plot(linspace(t1(r),t2(r),length(signal_normal)),uplim*amp_vec(k)+100*k,'Color',[0 0 0],'LineWidth',1)
                        plot(linspace(t1(r),t2(r),length(signal_normal)),lowlim*amp_vec(k)+100*k,'Color',[0 0 0],'LineWidth',1)

                        %                         
                        xlabel('Time')                    
                else
                        plot(linspace(t1(r),t2(r),length(signal_normal)),signal_normal*amp_vec(k)+100*k+500,'Color',colorvec(k,:))
                        hold on
                        plot(linspace(t1(r),t2(r),length(signal_normal)),uplim*amp_vec(k)+100*k+500,'Color',[0 0 0],'LineWidth',1)
                        plot(linspace(t1(r),t2(r),length(signal_normal)),lowlim*amp_vec(k)+100*k+500,'Color',[0 0 0],'LineWidth',1)

if (k==1 && lab==2)  % Threshold on PFC raw for most animals except rat12.
        %xo
    fact=2;    
    uplim=mean(signal_normal)+fact*std(signal_normal);
    lowlim=mean(signal_normal)-fact*std(signal_normal);
    %xo
    nsig=envelope(signal_normal);
    %nsig=(signal_normal>uplim(1)| signal_normal<lowlim(1));
    nsig=(nsig>uplim(1)| nsig<lowlim(1));
    nonZeroElements = nsig ~= 0;
    %minSeparation = 50000;
    %minSeparation = 10000;
    minSeparation=60*(fs_new); %10 seconds

    nonZeroElements = ~bwareaopen(~nonZeroElements, minSeparation);
    [labeledRegions, numRegions] = bwlabel(nonZeroElements);
    labeledRegions(labeledRegions~=0)=1;
    
    %Minimum of 3 sec duration.
    labeledRegions=min_duration(labeledRegions,3,fs_new);
    
    h=area(linspace(t1(r),t2(r),length(signal_normal)),labeledRegions.*1100);
    h.FaceColor=[1 0 0];
    h.FaceAlpha=0.2;
    % % alpha(0.3)

    % 
    % sh=stairs(labeledRegions);
    % 
    % bottom = min(sh.YData);
    % x = [sh.XData(1),repelem(sh.XData(2:end),2)];
    % y = [repelem(sh.YData(1:end-1),2),sh.YData(end)];
    % % plot(x,y,'y:') %should match stair plot
    % % Fill area
    % fill([x,fliplr(x)],[y,bottom*ones(size(y))], 'r')


    % Nsig=or(ischange(double(nsig)),nsig);

    %[nsig]=min_duration(nsig,1,fs_new);
    % bar(linspace(t1(r),t2(r),length(signal_normal)),nsig.*1100,'r'); alpha(0.6); 
    %xo
end

                        xlabel('Time')                    
                end
                clear uplim lowlim
            else
                %% Epoching data
                T1=1;
                T2=T1+2*(fs_new);
                con=1;
                while T2 <= length(signal_normal)
                %T2=T1+2*(300);    
                  new_seg = signal_normal(T1:T2);
                  NC(con,:)=new_seg;

                con=con+1;
                T2=T2+2*(fs_new);
                T1=T1+2*(fs_new);
                end
            [pxx,f]= periodogram(NC.',hann(size(NC.',1)),size(NC.',1),fs_new);    
            %%
            px=mean(pxx,2);
            s=semilogy(f,(px),'LineWidth',2);
            xlabel('Frequency (Hz)')
            ylabel('Power')
            xo
            printing(strcat(label1{lab},'_periodogram'))

            end
        % ylabel('Rat')
        end
        %xo
        end
            if plot_allanimals==1;
%                title(strcat('Normalized',{' '}, 'recordings',{' '},'Rat',num2str(rats(r))))
                title(['Rat',' ',num2str(rats(r)),' (',rat_condition{r},')'],'FontSize',12)

                yticks([100 200 300 400 500 600 700 800 900 1000])
%                 yticklabels({FN{1},FN{2},FN{3},FN{4},FN{5},FN{1},FN{2},FN{3},FN{4},FN{5}})
%yticklabels({['HPC',' ','Raw'],['HPC',' ','Delta'],['HPC',' ','Theta'],['HPC',' ','Spindles'],['HPC',' ','HighGamma'],['PFC',' ','Raw'],['PFC',' ','Delta'],['PFC',' ','Theta'],['PFC',' ','Spindles'],['PFC',' ','HighGamma']}) 
yticklabels({['HPC',' ','Raw'],['HPC',' ','Delta',' ','(0.1-4Hz)'],['HPC',' ','Theta',' ','(4-8Hz)'],['HPC',' ','Spindles',' ','(10-20Hz)'],['HPC',' ','HighGamma',' ','(80-250Hz)'],['PFC',' ','Raw'],['PFC',' ','Delta',' ','(0.1-4Hz)'],['PFC',' ','Theta',' ','(4-8Hz)'],['PFC',' ','Spindles',' ','(10-20Hz)'],['PFC',' ','HighGamma',' ','(80-250Hz)']}) 
                
                ylim([0 600+500])
                xlim([min(t1) max(t2)]);
                cd(selpath)
%                xo
%                 save(['stats_' 'Rat' num2str(rats(r)) '.mat'],'stats_vec','stats_vec_norm');
                %xo
                saveas(gcf,strcat('Rat',num2str(rats(r)),'_states60','.png'))
%                 savefig(gcf,strcat('xRat',num2str(rats(r)),'_all_bands_STD','.fig'),'compact') 
                close all
            end
    end
%end



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