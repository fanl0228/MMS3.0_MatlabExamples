

tic
clearvars
close all
setenv('CASCADE_SIGNAL_PROCESSING_CHAIN_MIMO', 'D:\ti\mmwave_studio_03_00_00_14\mmWaveStudio\MMS3.0_MatlabExamples\Hybrid_Beamforming')

UsedFrames = 10;
NumRx = 16;
num_chirps = 252;
LOG_ON = 0;
PLOT_ON =0;


dataFolder_Path = 'I:\20230918_EXP\TXBF_AngleSweep_Vib_Range2m_Angle+45_Sweep121\'; 

i = 1;
est_PSINR_all =[];
for Txbeam_angle = -60:60
    %Txbeam_angle = 60; % folderIdx = 121
    
    
    [est_PSINR, Txbeam_angle] = Each_Steering_Calculate_pSINR(dataFolder_Path, Txbeam_angle, UsedFrames, NumRx, num_chirps, LOG_ON, PLOT_ON);
    
    est_PSINR_all(i) = est_PSINR;
    
    i = i +1;
    toc
end

est_PSINR_all;

