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
selpath = uigetdir('C:\','Select folder with Downsampled_data');

plot_allanimals=0;
%Starting times for Rat2-Rat5.
t1(1) = duration([12 0 0]);
t1(2) = duration([11 53 0]);
t1(3) = duration([12 0 0]);
t1(4) = duration([11 54 0]);

t2=t1+seconds(durations_trials)'; %Values calculated from a previous iteration.
%% HPC
for r=1:length(rats)

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
load(cfold{1}) %HPC

HPC_normal=(HPC_downsampled-mean(HPC_downsampled))/std(HPC_downsampled);

if plot_allanimals==1;
    plot(linspace(t1(r),t2(r),length(HPC_normal)),HPC_normal+100*r)
    hold on
    xlabel('Time')
else
    %% Epoching data
    T1=1;
    T2=1;
    con=1;
    while T2 <= length(HPC_normal)
    T2=T1+2*(300);    
      new_seg = HPC_normal(T1:T2);
      NC(con,:)=new_seg;

    con=con+1;
    T1=T2;
    end
[pxx,f]= periodogram(NC.',hann(size(NC.',1)),size(NC.',1),300);    
%%
px=mean(pxx,2);
s=semilogy(f,(px),'LineWidth',2);
xlabel('Frequency (Hz)')
ylabel('Power')
printing('HPC_periodogram')
    
end
% ylabel('Rat')
end
    if plot_allanimals==1;
        title('Normalized HPC recording')
        yticks([100 200 300 400])
        yticklabels({'Rat 2','Rat 3','Rat 4','Rat 5'})
        printing('HPC')
    end
%% PFC

for r=1:length(rats)

cd('C:\Users\addri\Documents\Donders\Projects\CBD\Downsampled_data')

% %Find and access proper folder.
cfold=dir;
cfold={cfold.name};
cf1=cfold(cellfun(@(x) ~isempty(strfind(x,num2str(rats(r)))),cfold));
cd(cf1{1})

% %Find and load .mat files.
cfold=dir;
cfold={cfold.name};
cfold=cfold(cellfun(@(x) ~isempty(strfind(x,'.mat')),cfold));
load(cfold{2}) %PFC

PFC_normal=(PFC_downsampled-mean(PFC_downsampled))/std(PFC_downsampled);

    if plot_allanimals==1;
        plot(linspace(t1(r),t2(r),length(PFC_normal)),PFC_normal+100*r)
        hold on
        xlabel('Time')
        % ylabel('Rat')
    else
        %% Epoching data
        T1=1;
        T2=1;
        con=1;
        while T2 <= length(PFC_normal)
        T2=T1+2*(300);    
          new_seg = PFC_normal(T1:T2);
          NC(con,:)=new_seg;

        con=con+1;
        T1=T2;
        end
    [pxx,f]= periodogram(NC.',hann(size(NC.',1)),size(NC.',1),300);    
    %%
    px=mean(pxx,2);
    s=semilogy(f,(px),'LineWidth',2);
    xlabel('Frequency (Hz)')
    ylabel('Power')
    printing('PFC_periodogram')
    
end

    end
%end

    if plot_allanimals==1;
        title('Normalized PFC recording')
        yticks([100 200 300 400])
        yticklabels({'Rat 2','Rat 3','Rat 4','Rat 5'})
        printing('PFC')       
    end
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