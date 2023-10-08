

tic
clearvars
close all
setenv('CASCADE_SIGNAL_PROCESSING_CHAIN_MIMO', 'D:\ti\mmwave_studio_03_00_00_14\mmWaveStudio\MMS3.0_MatlabExamples\Hybrid_Beamforming')

UsedFrames = 10;
NumRx = 16;
num_chirps = 252;
LOG_ON = 1;
PLOT_ON =0;


dataFolder_Path = 'I:\20230928_Exp_NLOS\Indoor_Scene2\test1\'; 

i = 1;
est_PSINR_all =[];
Intensity_estSINR_all =[];
Phase_estSINR_all = [];
for Txbeam_angle = -60:60

    [Intensity_estSINR, Phase_estSINR, Txbeam_angle] = Each_Steering_Calculate_pSINR(dataFolder_Path, Txbeam_angle, UsedFrames, NumRx, num_chirps, LOG_ON, PLOT_ON);
    Intensity_estSINR_all(i) = Intensity_estSINR;
    Phase_estSINR_all(i) = Phase_estSINR;
    est_PSINR_all(i) = Intensity_estSINR + Phase_estSINR;
    
    i = i +1;
    toc
end

figure(220)
plot(Intensity_estSINR_all, 'r')
hold on;
plot(Phase_estSINR_all, 'b')
hold on;
plot(est_PSINR_all, 'g')
legend(["Intensity", "Phase", "All"]);

est_PSINR_all1=10*log10(10.^(Intensity_estSINR_all./10) + 10.^(Phase_estSINR_all./10));
hold on;
plot(est_PSINR_all1, 'k')
pause(0.01)
