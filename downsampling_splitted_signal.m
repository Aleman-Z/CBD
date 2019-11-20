function [PFC_downsampled]=downsampling_splitted_signal(LFP_PFC)

    LFP_PFC1=LFP_PFC(1:length(LFP_PFC)/2,1);
    LFP_PFC2=LFP_PFC(length(LFP_PFC)/2+1:end,1);

    clear LFP_PFC

    PFC_filt1=filtfilt(b,a,LFP_PFC1);
    clear LFP_PFC1
    PFC_downsampled1=downsample(PFC_filt1,fs/fs_new);
    clear PFC_filt1


    PFC_filt2=filtfilt(b,a,LFP_PFC2);
    clear LFP_PFC2
    PFC_downsampled2=downsample(PFC_filt2,fs/fs_new);
    clear PFC_filt2

    PFC_downsampled=[PFC_downsampled1; PFC_downsampled2];

end