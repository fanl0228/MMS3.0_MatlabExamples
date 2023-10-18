%%
% Binary Search 
%
% 采用二分查找的方式，搜索最优 TxBF 角度；
% 这里要求 RxBF = TxBF

tic
clearvars
close all
setenv('CASCADE_SIGNAL_PROCESSING_CHAIN_MIMO', 'D:\ti\mmwave_studio_03_00_00_14\mmWaveStudio\MMS3.0_MatlabExamples\Hybrid_Beamforming')

% Init Parameters of function Each_Steering_Calculate_pSINR 
dataFolder_Path = 'I:\LOS_Scene2_Raw\'; 
isStepAngle = 1;
TxBF_Angle = 0;
RxBF_Angle = TxBF_Angle;
UsedFrames = 20;
NumRx = 16;
num_chirps = 252;
Signal_FS = 2000;  % Doppler-Fs
LOG_ON = 1;
PLOT_ON = 0;
png_floder = strcat([dataFolder_Path(1:end-5), '_DRL\']);
if ~exist(png_floder,'dir')
    mkdir(png_floder)
end
LogFileId = fopen(strcat(png_floder, 'LogFile.txt'), 'w');


for SampleNum = 0:2
    for TxBF_Angle = -60: 60
   
        %SampleNum = floor(random('Uniform', 0, 3));
        %TxBF_Angle = floor(random('Uniform', -60, 61));
    
        FolderId = TxBF_Angle + 61;
        if (FolderId < 10)
            dataFolderName = strcat([dataFolder_Path, ...
                                    'BeamAngle00', num2str(FolderId), ...
                                    '_Sample0', num2str(SampleNum), ...
                                    '\']);
        else
            if (FolderId >= 10 && FolderId < 100)
                dataFolderName = strcat([dataFolder_Path, ...
                                        'BeamAngle0', num2str(FolderId), ...
                                        '_Sample0', num2str(SampleNum), ...
                                        '\']);
            else
                dataFolderName = strcat([dataFolder_Path, ...
                                        'BeamAngle', num2str(FolderId), ...
                                        '_Sample0', num2str(SampleNum), ...
                                        '\']);
            end
        end
        
        RxBF_Angle = TxBF_Angle;
        [Intensity_estSINR, Phase_estSINR, gframe_obj]...
                = Each_Steering_Calculate_pSINR(dataFolderName, ...
                                                isStepAngle, ...
                                                TxBF_Angle, ...
                                                RxBF_Angle, ...
                                                UsedFrames, ...
                                                NumRx, ...
                                                num_chirps, ...
                                                Signal_FS,...
                                                LOG_ON, ...
                                                PLOT_ON, ...
                                                LogFileId);
        
        Reward_PSINR = Intensity_estSINR + Phase_estSINR;
        disp([strcat('PSINR_Reward (dB):  ', num2str(Reward_PSINR))])
        
        %% -------------------------RL States Features
        % 1. Range-Doppler & Range-Angle, ObjRange-Angle
        for frameId = 1:length(gframe_obj)
            State_RD(frameId, :, :, :) = gframe_obj{1, frameId}.RangeDoppler_Profile;
            State_RA(frameId, :, :) = gframe_obj{1, frameId}.RangeAngle_Dynamic_Profile;
            State_RAObj(frameId, :) = squeeze(gframe_obj{1, frameId}.RA_Spectrum);
        end
        
        temp = split(dataFolderName,'\');
        angle_name = cell2mat(temp(end-1));
    
        State_RD_file = strcat(png_floder, angle_name, '_State_RD.mat');
        State_RA_file = strcat(png_floder, angle_name, '_State_RA.mat');
        State_RAObj_file = strcat(png_floder, angle_name, '_State_RAObj.mat');

        State_PSINR_file = strcat(png_floder, angle_name, '_Reward_PSINR.mat');
        State_Intensity_file = strcat(png_floder, angle_name, '_Reward_Intensity.mat');
        State_Phase_file = strcat(png_floder, angle_name, '_Reward_Phase.mat');

        save(State_RD_file, 'State_RD');
        save(State_RA_file, 'State_RA');
        save(State_RAObj_file, 'State_RAObj');

        save(State_PSINR_file, 'Reward_PSINR');
        save(State_Intensity_file, 'Intensity_estSINR');
        save(State_Phase_file, 'Phase_estSINR');

    
        % vis state information
        if PLOT_ON
            frameId_VIS = 1;
            figure(1)
            subplot(131)
            [y_axis, x_axis]=meshgrid(1: size(State_RD, 3), 1: size(State_RD, 2));
            fig2=surf(y_axis, x_axis, 10*log10(abs(squeeze(State_RD(frameId_VIS, :, :, 1)))));
            fig2.FaceColor = 'interp';
            fig2.LineStyle="none";
            colormap(gca,"jet");
            c = colorbar;
            c.Label.String = 'Relative Power(dB)';
            pause(0.05)
        
            subplot(132)
            [y_axis, x_axis]=meshgrid(1: size(State_RA, 3), 1: size(State_RA, 2));
            fig2=surf(y_axis, x_axis, 10*log10(abs(squeeze(State_RA(frameId_VIS, :, :)))));
            fig2.FaceColor = 'interp';
            fig2.LineStyle="none";
            colormap(gca,"jet");
            c = colorbar;
            c.Label.String = 'Relative Power(dB)';
            pause(0.05)
        
            subplot(133)
            plot(10*log10(abs(squeeze(State_RAObj(frameId_VIS, :)))));
            pause(0.05)
        end
        
        %TxBF_Angle = TxBF_Angle + 1;
        step_time = toc;
        disp(strcat('step processing time (s):  ', num2str(step_time)))
    end
end

fclose(LogFileId);
fclose all;


