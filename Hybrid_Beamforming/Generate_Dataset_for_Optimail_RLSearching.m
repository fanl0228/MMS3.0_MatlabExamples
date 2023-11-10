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

isStepAngle = 1;
TxBF_Angle = 0;
RxBF_Angle = TxBF_Angle;
UsedFrames = 10;
NumRx = 16;
num_chirps = 252;
Signal_FS = 2000;  % Doppler-Fs
LOG_ON = 0;
PLOT_ON = 0;

for LOSNUM = 20
    
    dataFolder_Path = strcat(['H:\NLOS', num2str(LOSNUM), '\']); 

    %dataFolder_Path = strcat(['I:\HyBF_Datasets\NLOS', num2str(LOSNUM), '\']); 
  

    png_floder = strcat([dataFolder_Path(1:end-1), '_DRL\']);
    if ~exist(png_floder,'dir')
        mkdir(png_floder)
    end
    LogFileId = fopen(strcat(png_floder, 'LogFile.txt'), 'w');
    
    if LOSNUM == 20
        SNUM = 2; 
    else
        SNUM = 100;
    end

    for SampleNum = SNUM
        
        rangeResolution = 0;   % calculate value
        est_PSINR_all =[];
        Intensity_estSINR_all =[];
        Phase_estSINR_all = [];
        range_angle_spec = [];
        step = 1;
        for TxBF_Angle = -60: 60
       
            % SampleNum = floor(random('Uniform', 0, 3));
            % TxBF_Angle = floor(random('Uniform', -60, 61));
        
            FolderId = TxBF_Angle + 61;
            if (FolderId < 10)
                if 1
                    dataFolderName = strcat([dataFolder_Path, ...
                                        'BeamAngle00', num2str(FolderId), ...
                                        '_Sample0', num2str(SampleNum), ...
                                        '\']);
                else
                    dataFolderName = strcat([dataFolder_Path, ...
                                        'Test_Vib_TXBF_BeamAngle00', num2str(FolderId), ...
                                        '\']);
                end
    
            else
                if (FolderId >= 10 && FolderId < 100)
                    if 1
                        dataFolderName = strcat([dataFolder_Path, ...
                                            'BeamAngle0', num2str(FolderId), ...
                                            '_Sample0', num2str(SampleNum), ...
                                            '\']);
                    else
                        dataFolderName = strcat([dataFolder_Path, ...
                                         'Test_Vib_TXBF_BeamAngle0', num2str(FolderId), ...
                                         '\']);
                    end
                else
                    if 1
                        dataFolderName = strcat([dataFolder_Path, ...
                                            'BeamAngle', num2str(FolderId), ...
                                            '_Sample0', num2str(SampleNum), ...
                                            '\']);
                    else
                        dataFolderName = strcat([dataFolder_Path, ...
                                         'Test_Vib_TXBF_BeamAngle', num2str(FolderId), ...
                                         '\']);
                    end
                end
            end
            
            RxBF_Angle = 0; %TxBF_Angle;
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
            
            % [Val_Intensity, idx_Intensity] = max(Intensity_estSINR);
            % Intensity_threshold = 0.7 * Val_Intensity;
            % for idx = 1:len(Intensity_estSINR)
            %     if Intensity_estSINR(idx) > Intensity_threshold
            % 
            %     end
            % end
    
            Reward_PSINR = Intensity_estSINR + Phase_estSINR;
            
            disp([strcat('PSINR_Reward (dB):  ', num2str(Reward_PSINR))])
            
            %% -------------------------RL States Features
            % 1. Range-Doppler & Range-Angle, ObjRange-Angle
            for frameId = 1:length(gframe_obj)
                State_RD_temp(frameId, :, :, :) = gframe_obj{1, frameId}.RangeDoppler_Profile;
                State_RA_temp(frameId, :, :) = gframe_obj{1, frameId}.RangeAngle_Dynamic_Profile;
                State_RAObj(frameId, :) = squeeze(gframe_obj{1, frameId}.RA_Spectrum);
    
    
                State_RDA_temp(frameId, :, :, :) = State_RD_temp(frameId, ...
                                                            1:floor(size(State_RD_temp,2)/2), ...   % range
                                                            [floor(size(State_RD_temp,3)/4)-1:floor(size(State_RD_temp,3)*3/4)], ... % doppler
                                                            :);
                State_RDA(frameId, :, :, :) = State_RDA_temp(frameId, :, :, :);
                % State_RDA(frameId, :, :, :) = fftshift(fft(State_RDA_temp(frameId, :, :, :), 64, 4), 4);
                % State_RDA(frameId, :, :, :) = single(10*log10(abs(State_RDA(frameId, :, :, :))));
    
    
                State_RA(frameId, :, :) = State_RA_temp(frameId, 1:floor(size(State_RD_temp,2)/2), :);
                
                State_RAObj(frameId, :) = State_RAObj(frameId, :);
    
            end
            
            % vis state information
            if PLOT_ON
                frameId_VIS = 1;
                figure(1)
                subplot(221)
                State_RD_avgAnt = squeeze(mean(State_RDA, 4));
                [y_axis, x_axis]=meshgrid(1: size(State_RD_avgAnt, 3), 1: size(State_RD_avgAnt, 2));
                fig2=surf(y_axis, x_axis, 10*log10(abs(squeeze(State_RD_avgAnt(frameId_VIS, :, :)))));
                %fig2=surf(y_axis, x_axis, squeeze(State_RD(frameId_VIS, :, :, 1)));
                fig2.FaceColor = 'interp';
                fig2.LineStyle="none";
                colormap(gca,"jet");
                c = colorbar;
                c.Label.String = 'Relative Power(dB)';
                title("Range-Doppler")
                pause(0.05)
            
                subplot(222)
                [y_axis, x_axis]=meshgrid(1: size(State_RA, 3), 1: size(State_RA, 2));
                fig2=surf(y_axis, x_axis, 10*log10(abs(squeeze(State_RA(frameId_VIS, :, :)))));
                fig2.FaceColor = 'interp';
                fig2.LineStyle="none";
                colormap(gca,"jet");
                c = colorbar;
                c.Label.String = 'Relative Power(dB)';
                title("Range-Angle")
                pause(0.05)
            
                subplot(223)
                plot(10*log10(abs(squeeze(State_RAObj(frameId_VIS, :)))));
                pause(0.05)
                
                State_RA_NoDopper = squeeze(mean(State_RDA, 3));
                subplot(224)
                [y_axis, x_axis]=meshgrid(1: size(State_RA_NoDopper, 3), 1: size(State_RA_NoDopper, 2));
                fig2=surf(y_axis, x_axis, 10*log10(abs(squeeze(State_RA_NoDopper(frameId_VIS, :, :)))));
                fig2.FaceColor = 'interp';
                fig2.LineStyle="none";
                colormap(gca,"jet");
                c = colorbar;
                c.Label.String = 'Relative Power(dB)';
                title("Range-Angle")
                pause(0.05)
    
            end % end if
    
            % Save data
            temp = split(dataFolderName,'\');
            angle_name = cell2mat(temp(end-1));
        
            State_RD_file = strcat(png_floder, angle_name, '_State_RDA.mat');
            State_RA_file = strcat(png_floder, angle_name, '_State_RA.mat');
            State_RAObj_file = strcat(png_floder, angle_name, '_State_RAObj.mat');
    
            State_PSINR_file = strcat(png_floder, angle_name, '_Reward_PSINR.mat');
            State_Intensity_file = strcat(png_floder, angle_name, '_Reward_Intensity.mat');
            State_Phase_file = strcat(png_floder, angle_name, '_Reward_Phase.mat');
    
            save(State_RD_file, 'State_RDA');
            save(State_RA_file, 'State_RA');
            save(State_RAObj_file, 'State_RAObj');
    
            save(State_PSINR_file, 'Reward_PSINR');
            save(State_Intensity_file, 'Intensity_estSINR');
            save(State_Phase_file, 'Phase_estSINR');
            
            %TxBF_Angle = TxBF_Angle + 1;
            step_time = toc;
            disp(strcat('step processing time (s):  ', num2str(step_time)))
        
        Intensity_estSINR_all(step) = Intensity_estSINR;
        Phase_estSINR_all(step) = Phase_estSINR;
        est_PSINR_all(step) = Intensity_estSINR + Phase_estSINR;
        
        if ~isempty(gframe_obj)
            ra_spec = abs(gframe_obj{1}.RA_Spectrum);
            for j = 2:length(gframe_obj)
                ra_spec = ra_spec + abs(gframe_obj{j}.RA_Spectrum);
            end
            ra_spec = ra_spec/length(gframe_obj);
            ra_spec_log = 10*log10(ra_spec);
            range_angle_spec(:, step) = ra_spec_log;
        else
            range_angle_spec(:, step) = zeros(1, 128); % range bin size = 128
        end
    
        rangeResolution = gframe_obj{1}.rangeResolution;
    
        step = step +1;
        
        
        end % end for TxBF_Angle = -60: 60

        %% -----------------------PSINR
        fig225=figure();
        search_angle = -floor(length(Intensity_estSINR_all)/2):1:floor(length(Intensity_estSINR_all)/2);
        %plot(search_angle, Intensity_estSINR_all, 'r')
        % hold on;
        plot(search_angle, Phase_estSINR_all, 'b')
        % hold on;
        % plot(search_angle, est_PSINR_all, 'g')
        hold off;
        % calculate the beamangle
        [val, idx] = max(est_PSINR_all);
        BeamAngle = idx - 61;
        %
        xlim([-floor(length(Intensity_estSINR_all)/2), floor(length(Intensity_estSINR_all)/2)])
        legend({"Intensity SINR", "Motion SINR", "PSINR"}, 'Location','best');
        xlabel("Beamsearching Angle")
        ylabel("SINR(dB)")
        xticks(-60:5:60)
        grid on;
        title({"Evaluate PSINR Metric", strcat("optimal BeamAngle: ", num2str(BeamAngle))})
        pause(0.01)
        % Save fig
        
        psinr_file = strcat([png_floder, 'LOS', num2str(LOSNUM),'_Sample', num2str(SampleNum),'_PSINR.png']);
        saveas(fig225, psinr_file, 'png');
        psinr_filefig = strcat([png_floder, 'LOS', num2str(LOSNUM),'_Sample', num2str(SampleNum),'_PSINR.fig']);
        saveas(fig225, psinr_filefig, 'fig');

        fig226=figure();
        sizeMarker = 50;
        colorMarker = linspace(1,10,length(Intensity_estSINR_all)); 
        scatter(Intensity_estSINR_all, Phase_estSINR_all, sizeMarker, colorMarker, 'filled');
        hold on;
        Intensity_max = max(Intensity_estSINR_all);
        Intensity_min = min(Intensity_estSINR_all);
        Phase_estSINR_max = max(Phase_estSINR_all);
        Phase_estSINR_min = min(Phase_estSINR_all);
        plot([Intensity_max*0.8, Intensity_max*0.8],[Phase_estSINR_min*0.8, Phase_estSINR_max*1.2], 'r--',LineWidth=1)
        xlim([Intensity_min*0.9, Intensity_max*1.1])
        ylim([Phase_estSINR_min*0.9, Phase_estSINR_max*1.1])
        xlabel("SINR_{cs} (dB)")
        ylabel("SINR_{bs} (dB)")
        
        fprintf(LogFileId, '================================================>>\n');
        fprintf(LogFileId, '=Optimail TxBeamforming Angle (deg): ====>>>>>>%s deg: \n', num2str(BeamAngle));

        
        plot(search_angle, 10.^(Phase_estSINR_all/10) / max(10.^(Phase_estSINR_all/10)), 'b')
    end

    
end
fclose(LogFileId);
fclose all;


