%  Copyright (C) 2020 Texas Instruments Incorporated - http://www.ti.com/
%
%
%   Redistribution and use in source and binary forms, with or without
%   modification, are permitted provided that the following conditions
%   are met:
%
%     Redistributions of source code must retain the above copyright
%     notice, this list of conditions and the following disclaimer.
%
%     Redistributions in binary form must reproduce the above copyright
%     notice, this list of conditions and the following disclaimer in the
%     documentation and/or other materials provided with the
%     distribution.
%
%     Neither the name of Texas Instruments Incorporated nor the names of
%     its contributors may be used to endorse or promote products derived
%     from this software without specific prior written permission.
%
%   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
%   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
%   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
%   A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
%   OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
%   SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
%   LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
%   DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
%   THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
%   (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
%   OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
%
%

tic
clearvars
close all
setenv('CASCADE_SIGNAL_PROCESSING_CHAIN_MIMO', 'D:\ti\mmwave_studio_03_00_00_14\mmWaveStudio\MMS3.0_MatlabExamples\Hybrid_Beamforming')


PLOT_ON = 1; % 1: turn plot on; 0: turn plot off
LOG_ON = 1; % 1: log10 scale; 0: linear scale
% numFrames_toRun = 10; %number of frame to run, can be less than the frame saved in the raw data
SAVEOUTPUT_ON = 0;
PARAM_FILE_GEN_ON = 1;
DISPLAY_RANGE_AZIMUTH_DYNAMIC_HEATMAP = 1 ; % Will make things slower
DEBUG_PLOTS = 0;                % optionally display debug plots while calibration data is being processed 


pro_path = getenv('CASCADE_SIGNAL_PROCESSING_CHAIN_MIMO');
input_path = strcat(pro_path,'\utils\cascade_json_parser\');
dataPlatform = 'TDA2'


% parameter file name for the test
pathGenParaFile = [input_path,'generateClibrationMatrix_param.m'];
% important to clear the same.m file, since Matlab does not clear cache
% automatically
clear(pathGenParaFile);

%module_param_file defines parameters to init each signal processing
%module 
% module_param_antennaCalib.m
module_param_file = [input_path, 'module_param.m']; 

searchBinsSkip = 64;            % number of bins to skip when looking for peak from corner reflector 
Start_Freq_GHz = 77;
TI_Cascade_Antenna_DesignFreq = 76.8; % antenna distance is designed for this frequency
Slope_MHzperus = 50;
targetRange = 2.1;  % estimated corner-reflector target range (in meters) 


% datacolor = colormap(jet(20));

% loop through each folder in the calibration directory create 
% matrix with dimensions [devices, number of TX, number of phase shift offsets] 
% with the phase values for the detected range-bin peak
% folder holding all of the datasets
dataFolder_Path = 'I:\20230918_EXP\TXBF_AngleSweep_Vib_Range4m_Angle+14.5_Sweep121\'; 
dataFolder = dir(dataFolder_Path);

%UsedFrame_Skip = 3; 

% find the first folder of cal data
for folderIdx = 1:length(dataFolder)
    if(dataFolder(folderIdx).isdir && dataFolder(folderIdx).name ~= "." && dataFolder(folderIdx).name ~= "..") 
        folderIdxStart = folderIdx;
        break;
    end

end
folderIdx = folderIdxStart; % start at first index that is a folder

%
%----------------------------------------get fundationmal informations of datasets to init. list------------------
% [beams, frames, samples, chirps, nRx]
if(dataFolder(folderIdx).isdir) 
    disp(['dataFolder(folderIdx).name = ', dataFolder(folderIdx).name]);

    % create folder structure for this iteration        
    fileNameCascade.dataFolderName = strcat([dataFolder(folderIdx).folder, '\', dataFolder(folderIdx).name, '\']);
    [fileIdx_unique] = getUniqueFileIdx(fileNameCascade.dataFolderName);
    
    [fileNameStruct] = getBinFileNames_withIdx(fileNameCascade.dataFolderName, fileIdx_unique{1});

    [numValidFrames, dataFileSize] = getValidNumFrames(fullfile(fileNameCascade.dataFolderName, fileNameStruct.masterIdxFile));

    %generate parameter file for the dataset
    parameter_file_gen_AngleSweep_json(fileNameCascade.dataFolderName, module_param_file, pathGenParaFile, dataPlatform);

    %generate calibration matrix object for the dataset
    genCalibrationMatrixObj = genCalibrationMatrixCascade('pfile', pathGenParaFile, ...
                                                          'calibrateFileName',fileNameCascade.dataFolderName, ...
                                                          'targetRange', targetRange);
    genCalibrationMatrixObj.binDataFile = fileNameStruct;
   
end


% simTopObj is used for top level parameter parsing and data loading and saving
simTopObj                   = simTopCascade('pfile', pathGenParaFile);
rangeFFTObj                 = rangeProcCascade('pfile', pathGenParaFile);
DopplerFFTObj               = DopplerProcClutterRemove('pfile', pathGenParaFile);
CFAR_DetectionObj           = CFAR_CASO('pfile', pathGenParaFile);
DOAObj                      = DOACascade('pfile', pathGenParaFile);

% get system level variables
platform                    = simTopObj.platform;
numValidFrames              = simTopObj.totNumFrames;
frameCountGlobal            = 0;
UsedFrames                  = 7;
if (UsedFrames > numValidFrames)
    UsedFrames              = numValidFrames;
end
dataFolder_Skip             = 1;


% 变量获取  
numTheta = ceil((length(dataFolder)-2)/dataFolder_Skip);
numFrames = UsedFrames; % ceil(numValidFrames/UsedFrame_Skip);
numSamples = genCalibrationMatrixObj.numSamplePerChirp;
numchirps = genCalibrationMatrixObj.nchirp_loops;
DopplerFFTSize = 2^(ceil(log2(numchirps)));
numRx = genCalibrationMatrixObj.numRxToEnable;

gBeam_Range_Profile = zeros(numTheta, numFrames, numSamples, numchirps, numRx);
gBeam_RangeDoppler_Profile = zeros(numTheta, numFrames, numSamples, DopplerFFTSize, numRx);

gBeams_Cube_Angle = zeros(numTheta, numFrames*numchirps, numRx);


% 定义每一帧数据所需变量
% detection_results_all_beams = {};
gBeam_frame_obj_rangeBins= [];
gBeam_dynamic_range_beam_spectrum = [];


% loop through all beam steering angles
for idxTheta = 1:ceil((length(dataFolder)-2)/dataFolder_Skip)
    if (folderIdx > length(dataFolder))
        fprintf(2, "folderIdx is out of the number of dataFolder_calib_data.. %d\n", folderIdx)
        break;
    end

    % 读取每个 TXBF 扫描角度的 数据文件夹
    if(dataFolder(folderIdx).isdir) 
        disp(['dataFolder(folderIdx).name = ', dataFolder(folderIdx).name]);

        % create folder structure for this iteration        
        fileNameCascade.dataFolderName = strcat([dataFolder(folderIdx).folder, '\', dataFolder(folderIdx).name, '\']);
        [fileIdx_unique] = getUniqueFileIdx(fileNameCascade.dataFolderName);
        
        [fileNameStruct] = getBinFileNames_withIdx(fileNameCascade.dataFolderName, fileIdx_unique{1});

        [numValidFrames, dataFileSize] = getValidNumFrames(fullfile(fileNameCascade.dataFolderName, fileNameStruct.masterIdxFile));

        %generate parameter file for the dataset
        parameter_file_gen_AngleSweep_json(fileNameCascade.dataFolderName, module_param_file, pathGenParaFile, dataPlatform);

        %generate calibration matrix object for the dataset
        genCalibrationMatrixObj = genCalibrationMatrixCascade('pfile', pathGenParaFile, ...
                                                              'calibrateFileName',fileNameCascade.dataFolderName, ...
                                                              'targetRange', targetRange);
        genCalibrationMatrixObj.binDataFile = fileNameStruct;
       
    end
    
    % 定义每一帧数据所需变量
    detection_results_all_frames = {};
    cnt = 1;
    gPhase_obj_validFrameNum = 0;

    % 读取每一帧数据
    for frameId = 1:UsedFrames
        % use second frame for calibration 
        genCalibrationMatrixObj.frameIdx = frameId; 
        frameCountGlobal = frameCountGlobal+1;
        

        % force cal object to single TX phase/frame (numChirpPerLoop <-> numTXPhases)
        %genCalibrationMatrixObj.nchirp_loops = 1;
        genCalibrationMatrixObj.TxToEnable = 1;
        
        % read in .bin files， rawData = [samples, chirps, Rx]
        rawData = cascade_Read_rawData(genCalibrationMatrixObj);
        % if mod(frameId, 10)==1
        %     fprintf('Processing %3d frame...\n', frameId)
        % end
        
        %% -----------------Range Doppler FFT------------------------------
        %perform 2D FFT
        rangeFFTOut = [];
        DopplerFFTOut = [];
        
        % range FFT
        rangeFFTOut(:,:,:)     = datapath(rangeFFTObj, rawData(:,:,:));
        % visualization Range FFT Output
        if(0)
            figure(2)
            plot(10*log10(abs(rangeFFTOut(:,:,1))));
            title("Plot Range FFT Output");
            pause(0.01);
        end
        gBeam_Range_Profile(idxTheta, frameId,:,:,:) = rangeFFTOut(:,:,:);
        
        % Doppler FFT
        DopplerFFTOut(:,:,:)   = datapath(DopplerFFTObj, rangeFFTOut(:,:,:));
        gBeam_RangeDoppler_Profile(idxTheta, frameId,:,:,:) = DopplerFFTOut(:,:,:);
        % Visualization Doppler FFT Output
        if(0)
            figure(2)
            fig2 = pcolor(10*log10( abs( DopplerFFTOut(:,:,1) ) ) );
            fig2.FaceColor = 'interp';
            fig2.LineStyle="none";
            title("Imagesc Range-Doppler FFT Output");
            pause(0.01)
        end
        
        sig_integrate = 10*log10(sum((abs(DopplerFFTOut)).^2, 3) + 1);
        
        % CFAR Detection
        detection_results = datapath(CFAR_DetectionObj, DopplerFFTOut);
        % detection_results_all_frames{cnt} =  detection_results;  % replaced by angleEst 

        detect_all_points = [];
        for iobj = 1:length(detection_results)
            detect_all_points(iobj, 1)=detection_results(iobj).rangeInd;
            detect_all_points(iobj, 2)=detection_results(iobj).dopplerInd_org;
            detect_all_points(iobj, 4)=detection_results(iobj).estSNR;
        end
        
       

        %% Range - angle, point cloud generation.
        angles_all_points = [];
        xyz = [];
        if ~isempty(detection_results)
            % DOA, the results include detection results + angle estimation results.
            % access data with angleEst{frame}(objectIdx).fieldName
            Txbeam_angle = (idxTheta-1)*dataFolder_Skip - floor((length(dataFolder)-2)/2);
            angleEst = datapath(DOAObj, detection_results, Txbeam_angle, DopplerFFTOut);
             
            % Visualization
            if PLOT_ON
                figure(3);
                set(gcf,'units','normalized','outerposition',[0.1 0.1 0.8 0.8])                
                
                subplot(2,3,1)               
                plot((1:size(sig_integrate,1))*CFAR_DetectionObj.rangeBinSize, sig_integrate(:,size(sig_integrate,2)/2+1),'g','LineWidth',4);
                hold on; grid on;
                for ii=1:size(sig_integrate,2)
                    
                    if ~isempty(detection_results)
                        ind = find(detect_all_points(:,2) == ii);
                        if (~isempty(ind))
                            plot((1:size(sig_integrate,1))*CFAR_DetectionObj.rangeBinSize, sig_integrate(:, ii));
                            hold on; grid on;
    
                            rangeInd = detect_all_points(ind,1);
                            plot(rangeInd*CFAR_DetectionObj.rangeBinSize, sig_integrate(rangeInd, ii), ...
                                    'o','LineWidth',2,...
                                    'MarkerEdgeColor','k',...
                                    'MarkerFaceColor',[.49 1 .63],...
                                    'MarkerSize',10);
                        end
                    end
                end
                xlabel('Range(m)');
                ylabel('Receive Power (dB)')
                title(['Range Profile(zero Doppler - thick green line): frameID ' num2str(frameId)]);
                hold off;
   
                subplot(2,3,2);
                %subplot_tight(2,2,2,0.1)
                sig_integrate_noDC = sig_integrate(:, [1:size(sig_integrate,2)/2-1 size(sig_integrate,2)/2+3:end]);
                % [Y,X]=meshgrid(1:1:size(sig_integrate_noDC,1), -size(sig_integrate_noDC,2)/2+1:1:size(sig_integrate_noDC,2)/2);
                % Y = Y * CFAR_DetectionObj.rangeBinSize;
                % X = X * CFAR_DetectionObj.velocityBinSize;
                fig2 = pcolor(sig_integrate_noDC);
                fig2.FaceColor = 'interp';
                fig2.LineStyle="none";
                colormap(gca,"jet");
                c = colorbar;
                c.Label.String = 'Relative Power(dB)';
    
                [RD_SNR_max_val, RD_SNR_max_idx]= max(detect_all_points(:, 4));
                title({' Range/Velocity Plot No-DC', ...
                        strcat('obj1 SNR=', num2str(detection_results(1).estSNR)) ...
                        strcat(['SNR_{max}(dB)' num2str(RD_SNR_max_val) ...
                                ' Range(m)=' num2str(detect_all_points(RD_SNR_max_idx, 1)*CFAR_DetectionObj.rangeBinSize)] )  ...
                        });
                pause(0.01)
            
            end

            if ~isempty(angleEst)
                %% Filting using the doppler, estSNR, and Range angle fft
                % 去掉doppler_corr为0的目标,并且SNR取 top10
                idxs1 = find([angleEst.doppler_corr]);
                
                [~, idxs2] = sort([angleEst.estSNR], "descend");
                if (length(idxs2) > 10)
                    SNR_TOPN = 10;
                else
                    SNR_TOPN = length(idxs2);
                end
                idxs3 = intersect(idxs1, idxs2);

                % 
                % 计算target信噪比，从中取SNR 前top N个点作为基础点。 angleEst
                for iobj = 1:length(idxs3)
                    out_frame(1, iobj).rangeInd                    = angleEst(1, iobj).rangeInd;
                    out_frame(1, iobj).dopplerInd                  = angleEst(1, iobj).dopplerInd;
                    out_frame(1, iobj).range                       = angleEst(1, iobj).range;
                    out_frame(1, iobj).doppler_corr                = angleEst(1, iobj).doppler_corr;
                    out_frame(1, iobj).dopplerInd_org              = angleEst(1, iobj).dopplerInd_org;
                    out_frame(1, iobj).noise_var                   = angleEst(1, iobj).noise_var;
                    out_frame(1, iobj).bin_val                     = angleEst(1, iobj).bin_val;
                    out_frame(1, iobj).estSNR                      = angleEst(1, iobj).estSNR;
                    out_frame(1, iobj).doppler_corr_overlap        = angleEst(1, iobj).doppler_corr_overlap;
                    out_frame(1, iobj).doppler_corr_FFT            = angleEst(1, iobj).doppler_corr_FFT;
                    out_frame(1, iobj).angles                      = angleEst(1, iobj).angles;
                    out_frame(1, iobj).spectrum                    = angleEst(1, iobj).spectrum;
                    out_frame(1, iobj).range_beam_spectrum         = angleEst(1, iobj).range_beam_spectrum;
                end
                detection_results_all_frames{cnt} = out_frame;
                
                % calculate target azimuth
                [~, idx_angle] = max( 10*log10( abs( [angleEst.range_beam_spectrum] ) ) );
                if(ismember(idx_angle, idxs3))
                    obj_rangebin = idx_angle;
                else
                    obj_rangebin= [detection_results_all_frames{cnt}.rangeInd];
                end

                % [Beam, frame, objs_range_bin]
                gBeam_frame_obj_rangeBins(idxTheta, frameId, 1:length(obj_rangebin)) = obj_rangebin; 
                gBeam_frame_obj_estSNR(idxTheta, frameId, 1:length(obj_rangebin))  = detection_results_all_frames{cnt}.estSNR;

                gBeam_dynamic_range_beam_spectrum(idxTheta, frameId, :)  = detection_results_all_frames{cnt}.range_beam_spectrum;
                
                gPhase_obj_validFrameNum = gPhase_obj_validFrameNum + 1;


                % plot 3D point cloud
                if length(out_frame) > 0
                    for iobj = 1:length(out_frame)
                        angles_all_points(iobj,1:2)=out_frame(iobj).angles(1:2);
                        angles_all_points(iobj,3)=out_frame(iobj).estSNR;
                        angles_all_points(iobj,4)=out_frame(iobj).rangeInd;
                        angles_all_points(iobj,5)=out_frame(iobj).doppler_corr;
                        angles_all_points(iobj,6)=out_frame(iobj).range;
                        %switch left and right, the azimuth angle is flipped
                        xyz(iobj,1) = angles_all_points(iobj,6) * sind( angles_all_points(iobj,1) * -1 ) * cosd( angles_all_points(iobj,2) );
                        xyz(iobj,2) = angles_all_points(iobj,6) * cosd( angles_all_points(iobj,1) * -1 ) * cosd( angles_all_points(iobj,2) );
                        %switch upside and down, the elevation angle is flipped
                        xyz(iobj,3) = angles_all_points(iobj,6) * sind(angles_all_points(iobj,2) * -1);
                        xyz(iobj,4) = out_frame(iobj).doppler_corr;
                        xyz(iobj,9) = out_frame(iobj).dopplerInd_org;
                        xyz(iobj,5) = out_frame(iobj).range;
                        xyz(iobj,6) = out_frame(iobj).estSNR;
                        xyz(iobj,7) = out_frame(iobj).doppler_corr_overlap;
                        xyz(iobj,8) = out_frame(iobj).doppler_corr_FFT; 
                    end

                    angles_all_all{cnt} = angles_all_points;
                    xyz_all{cnt}  = xyz;
                    maxRangeShow = CFAR_DetectionObj.rangeBinSize * rangeFFTObj.rangeFFTSize;
                    %tic
                    
                    if PLOT_ON
                        moveID = find(abs(xyz(:,4))>=0);
                        subplot(2,3,3);                        
                        
                        if cnt==1
                            scatter3(xyz(moveID,1),xyz(moveID,2),xyz(moveID,3),45,(xyz(moveID,4)), 'filled');
                            
                        else
                            yz = [xyz_all{cnt}; xyz_all{cnt-1}];
                            scatter3(xyz(moveID,1),xyz(moveID,2),xyz(moveID,3),45,(xyz(moveID,4)),'filled');
                        end
                        colormap(gca,"jet");
                        c = colorbar;
                        c.Label.String = 'velocity (m/s)';                        
                        grid on;
                        
                        xlim([-20 20])
                        ylim([1 maxRangeShow])
                        %zlim([-4 4])
                        zlim([-5 5])
                        xlabel('X (m)')
                        ylabel('y (m)')
                        zlabel('Z (m)')                        
                        
                        view(2)                        
                        title(' 3D point cloud');
                        
                        % plot range and azimuth heatmap
                        subplot(2,3,4)
                        % STATIC_ONLY: 1 = plot heatmap for zero-Doppler; 0 = plot heatmap for nonzero-Doppler
                        STATIC_ONLY = 1; 
                        minRangeBinKeep = 2;
                        rightRangeBinDiscard =  20;
                        [mag_data_static(:,:,frameCountGlobal), ...
                        mag_data_dynamic(:,:,frameCountGlobal), ...
                        y_axis, x_axis] = plot_range_azimuth_2D(CFAR_DetectionObj.rangeBinSize, ...
                                                                DopplerFFTOut,...
                                                                length(genCalibrationMatrixObj.TxForMIMOProcess), ...
                                                                length(genCalibrationMatrixObj.RxForMIMOProcess), ...
                                                                CFAR_DetectionObj.antenna_azimuthonly, ...
                                                                LOG_ON, ...
                                                                STATIC_ONLY, ...
                                                                PLOT_ON, ...
                                                                minRangeBinKeep, ...
                                                                rightRangeBinDiscard);
                        title('range/azimuth heat map static objects')

                        if (DISPLAY_RANGE_AZIMUTH_DYNAMIC_HEATMAP)                   
                            subplot(2,3,5);
                            surf(y_axis, x_axis, (mag_data_static(:,:,frameCountGlobal)).^0.4,'EdgeColor','none');
                            colormap(gca,"jet");
                            view(2);
                            xlabel('meters');    
                            ylabel('meters');
                            title({'Static Range-Azimuth Heatmap',strcat('Current Total Frame Number = ', num2str(frameCountGlobal))})
                            
                            subplot(2,3,6);
                            surf(y_axis, x_axis, (mag_data_dynamic(:,:,frameCountGlobal)).^0.4,'EdgeColor','none');
                            colormap(gca,"jet");
                            view(2);    
                            xlabel('meters');    
                            ylabel('meters');
                            title({'Dynamic HeatMap', strcat('curent Sweep angle= ',num2str(Txbeam_angle))});
                            
                        end
                        pause(0.1) 

                    end % PLOT_ON
        
                end % if length(angleEst) > 0
            
                cnt = cnt + 1;  
            end % ~isempty(angleEst) 

        end % ~isempty(detection_results)
        
        %% 

    end % frameId


    %% connect all chirps, and extract phase processing 
    gPhase_RXs = [];
    for frameId = 1:gPhase_obj_validFrameNum %ceil(numValidFrames/UsedFrame_Skip)
        targetbin = gBeam_frame_obj_rangeBins(idxTheta, frameId, 1);
        if (targetbin == 0)
            continue;
        else
            gPhase_RXs = cat(1, gPhase_RXs, ...
                        squeeze(gBeam_Range_Profile(idxTheta, ...
                                                    frameId, ...
                                                    gBeam_frame_obj_rangeBins(idxTheta, frameId, 1), ...
                                                    :, :) ...
                                                    ));
            %gPhase_obj_validFrameNum = gPhase_obj_validFrameNum+1;
        end
    end
    
    %
    % ---------------------Processing All Beam Space one target phase ---------------------
    
    % IQ plan analysis
    if (0)
        figure(6)
        set(gcf,'units','normalized','outerposition',[0.1 0.1 0.8 0.8])
        sizeMarker = 10;
        colorMarker = colormap(parula(size(gPhase_RXs, 1)));
        fig6 = tight_subplot(2,8,[.05 .05],[.1 .05],[.05 .05]);
        for i=1:size(gPhase_RXs, 2)
            axes(fig6(i));
            scatter(real(gPhase_RXs(:, i)), imag(gPhase_RXs(:, i)), sizeMarker, colorMarker);
            colorbar;
            title(sprintf("Rx: %d", i));
        end        
        pause(0.01)
    end
    
    if (~isempty(gPhase_RXs))


        for idxRX = 1:genCalibrationMatrixObj.numRxToEnable
            % 1. get angle 
            Phase_Angle = unwrap(angle(gPhase_RXs(:, idxRX)));

            if (DEBUG_PLOTS)    
                figure(7);
                plot(Phase_Angle);
                title("unwrap angle value")
                pause(0.01)
                
                signal = diff(Phase_Angle);
                signal = signal(1:length(signal)-1);

                figure(8);
                plot(signal);
                title("[Comparing] original phase diff signal,, No Smooth")
                pause(0.01)
            end
            
            % 2.fram 之间的 gap 相位不连续优化处理 
            if (gPhase_obj_validFrameNum > 2)
                % step1 去除 由 number chirp 引起的 相位突变
                angle_gap = 0.0;
                for idxChirp = 1:size(Phase_Angle, 1)/genCalibrationMatrixObj.nchirp_loops-1 

                    angle_gap = angle_gap ...
                                + Phase_Angle(idxChirp * genCalibrationMatrixObj.nchirp_loops + 1) ...
                                - Phase_Angle(idxChirp * genCalibrationMatrixObj.nchirp_loops);
                    
                    Phase_Angle(idxChirp * genCalibrationMatrixObj.nchirp_loops : ...
                                (idxChirp+1) * genCalibrationMatrixObj.nchirp_loops - 1) ...
                                    = Phase_Angle(idxChirp * genCalibrationMatrixObj.nchirp_loops : ...
                                                (idxChirp + 1) * genCalibrationMatrixObj.nchirp_loops - 1) - angle_gap;
                end
                

                % step2
                for idxChirp = 1:size(Phase_Angle, 1)/genCalibrationMatrixObj.nchirp_loops-1 
                    th1 = abs( Phase_Angle(idxChirp*genCalibrationMatrixObj.nchirp_loops-1) ...
                                -  Phase_Angle(idxChirp*genCalibrationMatrixObj.nchirp_loops));
                    th2 = abs( Phase_Angle(idxChirp*genCalibrationMatrixObj.nchirp_loops+1) ...
                                -  Phase_Angle(idxChirp*genCalibrationMatrixObj.nchirp_loops));
                    
                    Threshold = abs(Phase_Angle(idxChirp*genCalibrationMatrixObj.nchirp_loops+1) ...
                                    - Phase_Angle(idxChirp*genCalibrationMatrixObj.nchirp_loops-1));
            
                    if th1 + th2 > Threshold
                        Phase_Angle(idxChirp*genCalibrationMatrixObj.nchirp_loops) ...
                                = 0.5*(Phase_Angle(idxChirp*genCalibrationMatrixObj.nchirp_loops+1) ...
                                    + Phase_Angle(idxChirp*genCalibrationMatrixObj.nchirp_loops-1));
                    end
                end
            
            end
            
            % 保存每一个RX 通道的相位数据，【考虑后续在RX通道做beamforming】
            gBeams_Cube_Angle(idxTheta, 1:length(Phase_Angle), idxRX) = Phase_Angle;
            
            if(0)
                % step 相位差分
                Beam_Sgnal = diff(Phase_Angle);
                Beam_Sgnal = Beam_Sgnal(1:length(Beam_Sgnal) - 1);
                figure(9)
                subplot(121)
                plot(unwrap(angle(gPhase_RXs(:, idxRX))),'color', 'blue' );
                hold on;
                plot(Phase_Angle,'color', 'red');
                hold off;
                title({"Phase Signal", strcat("Target Bin",num2str(targetbin))})
                subplot(122)
                plot(Beam_Sgnal);
                title("gBeams Phase Angle Diff")
                pause(0.01)
            end %if(DEBUG_PLOTS)
                
        end % idxRx    
    end %if (~isempty(gPhase_RXs))

    % 删除局部变量
    clearvars Cube_Angle;
    clearvars Beam_Sgnal;
    clearvars Beam_Sgnal1;
    clearvars Phase_Angle;
    
    % 可视化 所有 RX 的差分信号
    if (0)
        figure(10)
        fig10 = tight_subplot(2,8,[.05 .05],[.1 .05],[.05 .05]);
        for i=1:size(gBeams_Cube_Angle, 3)
            axes(fig10(i));
            plot(squeeze(gBeams_Cube_Angle(idxTheta, :, i)));
            title(["Rx Number: " num2str(i)]);
        end
    end
    
    %% ------------------------------------------对每个角度设置计算其 pSINR

    % 所有RX 天线求平均
    gBeams_Angle = mean(gBeams_Cube_Angle(idxTheta, :, :), 3);
    gBeams_Angle_diff = diff(gBeams_Angle);
    if (1)
        figure(11)
        subplot(121)
        plot(gBeams_Angle,'color', 'blue');
        title({"Phase Signal", strcat("Target Bin",num2str(targetbin))})
        subplot(122)
        plot(gBeams_Angle_diff);
        title("gBeams Phase Angle Diff")
        pause(0.01)
    end
    
    % Calculate Sensing Signal Intensity SNR
    % SNR part 1
    target = mode(gBeam_frame_obj_rangeBins(idxTheta, :, 1), 'all');
    
    target_idx = find(gBeam_frame_obj_rangeBins(idxTheta, :, 1) == target);

    Intensity_SNR(idxTheta) = sum(gBeam_frame_obj_estSNR(idxTheta, target_idx, 1), 2) / length(target_idx);
    
    disp(strcat("======>>>>> Target Range bin = ", num2str(target)))
    disp(strcat("======>>>>> Target Index Range bin = ", num2str(target_idx)))
    disp(strcat("======>>>>> Intensity_SNR(dB) = ", num2str(Intensity_SNR(idxTheta))))

    PhaseSingal_SNR = [];  




    folderIdx = folderIdx + dataFolder_Skip;
    toc
end

% Test RXBF result
if (1)
    spectrogram_AOD = squeeze( gBeam_dynamic_range_beam_spectrum(:, 2, :) );
    SweepAngles = (1:size(spectrogram_AOD, 1) * dataFolder_Skip) - ceil(size(spectrogram_AOD, 1)/2);

    sine_theta = sind(SweepAngles);
    cos_theta = sqrt(1-sine_theta.^2);
    indices_1D = 1:size(spectrogram_AOD, 2);

    indices_1D_range = (1:size(spectrogram_AOD, 2))*genCalibrationMatrixObj.rangeResolution;
    [R_mat, sine_theta_mat] = meshgrid(indices_1D_range, sine_theta);
    [~, cos_theta_mat] = meshgrid(indices_1D_range, cos_theta);
    x_axis = R_mat.*cos_theta_mat;
    y_axis = R_mat.*sine_theta_mat;

    figure(19)
    surf(y_axis, x_axis, abs(spectrogram_AOD).^0.2,'EdgeColor','none');
    %xlim([-5 5])
    %ylim([0 10]);
    view(2);
    xlabel('Meters')
    ylabel('Meters')

    % calculate target azimuth
    [val, idx] = max( 10*log10( abs( sum(spectrogram_AOD,2) ) ) );
    target_angle = idx * dataFolder_Skip - ceil(size(spectrogram_AOD, 1)/2);

    title({'RxBF stich range/azimuth', strcat("target angle (deg) = ", num2str(target_angle))})


end










% 删除原始数据文件相关的 不用局部变量
clearvars fileIdx_unique;
clearvars fileNameCascade;
clearvars fileNameStruct;
clearvars rawData;


%% 波束增益绘制
gBeamGain_power = zeros(1, numTheta);
figure(20)
subplot(1,2,1)
obj1_rangebin = gBeam_frame_obj_rangeBins(idxTheta, frameId, 1);
gBeamGain_power = gBeam_Range_Profile(:, 2, obj1_rangebin, size(gBeam_Range_Profile, 4)/2+1, 1);
gBeamGain_power = 10*log10(abs(gBeamGain_power));
plot(fliplr(gBeamGain_power));
title("gBeamGain Power Pattern");

subplot(1,2,2)
theta = ((1:dataFolder_Skip:(length(dataFolder)-2))- floor((length(dataFolder)-2)/2) + 90) /180*pi; % 极坐标
polarplot(theta, fliplr(gBeamGain_power/max(gBeamGain_power)), 'r-','linewidth',1);
rlim([0.8, 1])
pause(0.01)

% 删除处理相关的使用变量 gBeamGain_power
clearvars gBeamGain_power;


%% Range-Angle Stich
PLOT_DYNAMIC = 1;
ratio = 0.5;
if PLOT_DYNAMIC 
    DopplerPower = squeeze( mean( sum( mean( abs(gBeam_RangeDoppler_Profile(:, 2, :, :, :)), 5), 3), 1) );
    DopplerPower_noDC = DopplerPower([1: numchirps/2-1 numchirps/2+3:end]);
    [peakVal peakInd] = max(DopplerPower_noDC);
    threshold = peakVal*ratio;
    indSel = find(DopplerPower_noDC >threshold);
    for ii = 1:length(indSel)
        if indSel(ii) > DopplerFFTSize/2-1
            indSel(ii) = indSel(ii) + 3;
        end
    end
    angle_range_dynamic = squeeze( gBeam_RangeDoppler_Profile(:, 2, :, indSel, :) );
    gBeam_RangeDoppler = angle_range_dynamic;

else
    gBeam_RangeDoppler_Profile_zeroDop = squeeze(gBeam_RangeDoppler_Profile(:,2,:,numchirps/2+1,:));
    gBeam_RangeDoppler = gBeam_RangeDoppler_Profile_zeroDop;
end

SweepAngles=(1:dataFolder_Skip:(length(dataFolder)-2)) - ceil((length(dataFolder)-2)/2);
range_angle_stich = zeros(size(gBeam_RangeDoppler, 2), length(SweepAngles));

D_RX = DOAObj.D_RX;
dlambda = DOAObj.antDis;

for idxTheta = 1:length(SweepAngles)-1
    angleTX = SweepAngles(idxTheta);
    wx = sind(angleTX);
    a1_az = exp(1j*2*pi*dlambda*(D_RX*wx));
    
    for irange = 1:size(gBeam_RangeDoppler, 2)
        
        RX_data = squeeze(gBeam_RangeDoppler(idxTheta, irange, :));
        
        range_angle_stich(irange, idxTheta) = a1_az*(RX_data'*RX_data)*a1_az';
    end
end
sine_theta = sind(SweepAngles);
cos_theta = sqrt(1-sine_theta.^2);
indices_1D_half = (1:size(gBeam_RangeDoppler, 2)/2);

indices_1D_half_range = (1:size(gBeam_RangeDoppler, 2)/2)*genCalibrationMatrixObj.rangeResolution;
[R_mat, sine_theta_mat] = meshgrid(indices_1D_half_range, sine_theta);
[~, cos_theta_mat] = meshgrid(indices_1D_half_range, cos_theta);
x_axis = R_mat.*cos_theta_mat;
y_axis = R_mat.*sine_theta_mat;

range_angle_stich_half = (range_angle_stich(indices_1D_half + size(gBeam_RangeDoppler, 2)/2, :).');
figure(12)
surf(y_axis, x_axis, abs(range_angle_stich_half).^0.2,'EdgeColor','none');
%xlim([-5 5])
%ylim([0 10]);
view(2);
xlabel('Angle')
ylabel('Meters')
title('stich range/azimuth')

% --------- 删除变量
clearvars gBeam_RangeDoppler_Profile_zeroDop;
%clearvars gBeam_RangeDoppler_Profile;
clearvars range_angle_stich;
clearvars SweepAngles;
clearvars range_angle_stich_half;


%%
if (0)
    TargetBin_Vis = 10;
    Frame_Vis = 2;
    
    for i = 1:size(gBeam_RangeDoppler_Profile,1)
        figure(13+i)
        imagesc(10*log2(abs(squeeze(gBeam_RangeDoppler_Profile(i, Frame_Vis, :, :, 1)))));
        title("After 2D Range-Doppler FFT");
    end
end

TargetBin_Vis = 10;
Frame_Vis = 2;
sumDoppler(:,:) = 10*log2(abs(squeeze(sum(gBeam_RangeDoppler_Profile(:, Frame_Vis, :, :, 1),4))));
CM13= colormap(jet(size(sumDoppler,1)));
figure(111)
for i =1:size(sumDoppler,1)
    plot(sumDoppler(i,:), 'Color', CM13(i,:),LineWidth=1);
    hold on;
    leg_str{i} = ['Theta ', num2str(i)];
end
legend(leg_str);

%%  在 Beam 所扫描的平面上寻找最优区域，计算 PSINR

% step1. 计算beam噪声平面,在RX1
nFFT = 256;
DopplerFs = 1000;
[nBeams, nPhases, nRx] = size(gBeams_Cube_Angle);
phase_diff = zeros(nBeams, nPhases-1);
phase_PSD_Half = zeros(nBeams, nFFT/2);
wind = hann(nFFT);
for idxTheta = 1:nBeams
    phase_diff(idxTheta, :) = diff(gBeams_Cube_Angle(idxTheta, :, 1));

    phase_PSD = (abs(fftshift(fft(phase_diff(idxTheta, 1:nFFT).*wind', nFFT)))).^2/(DopplerFs*nFFT);

    phase_PSD_Half(idxTheta, :) = phase_PSD(length(phase_PSD)/2+1:length(phase_PSD));
end
% 计算相邻 beam 信号能量差
for idx = 1:nBeams-1
    phase_PSD_Half_diff(idx) = sum(phase_PSD_Half(idx+1, :)) / sum(phase_PSD_Half(idx, :));
end

[val, idx] = max(phase_PSD_Half_diff);

beam_angle= idx * dataFolder_Skip - ceil((length(dataFolder)-2)/2);

sprintf("The maximum value of adjacent beams angle: %d, value: %.2f", beam_angle, val)
%% 结果算的不准










toc
timeProcessing = toc;


