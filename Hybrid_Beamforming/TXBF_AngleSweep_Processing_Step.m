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
module_param_file = [input_path, 'module_param_antennaCalib.m']; 

searchBinsSkip = 64;            % number of bins to skip when looking for peak from corner reflector 
Start_Freq_GHz = 77;
TI_Cascade_Antenna_DesignFreq = 76.8; % antenna distance is designed for this frequency
Slope_MHzperus = 50;
targetRange = 2.1;  % estimated corner-reflector target range (in meters) 

DEBUG_PLOTS = 0;                % optionally display debug plots while calibration data is being processed 

if(DEBUG_PLOTS) 
    fig1 = figure(1); % Range FFT 
    axes1 = axes;
    
    fig2 = figure(2); % Range-Doppler FFT
    axes2 = axes;
    
    fig3 = figure(3); % Range-Angle FFT 
    axes3 = axes;
    
    fig4 = figure(4); % FFT phase (deg)
    axes4 = axes;

    fig5 = figure(5); % FFT magnitude (dB)
    axes5 = axes;
end
datacolor = colormap(jet(20));

% loop through each folder in the calibration directory create 
% matrix with dimensions [devices, number of TX, number of phase shift offsets] 
% with the phase values for the detected range-bin peak
% folder holding all of the datasets
dataFolder_Path = 'I:\20230918_EXP\TXBF_AngleSweep_Vib_Range2m_Angle+25.8_Sweep121\'; 
dataFolder_Skip = 4;
dataFolder = dir(dataFolder_Path);

UsedFrame_Skip = 3;

% find the first folder of cal data
for folderIdx = 1:length(dataFolder)
    if(dataFolder(folderIdx).isdir && dataFolder(folderIdx).name ~= "." && dataFolder(folderIdx).name ~= "..") 
        folderIdxStart = folderIdx;
        break;
    end

end
folderIdx = folderIdxStart; % start at first index that is a folder

%
%----------------------------------------get fundationmal informations of datasets------------------
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

numTheta = ceil((length(dataFolder)-2)/dataFolder_Skip);
numFrames = ceil(numValidFrames/UsedFrame_Skip);
numSamples = genCalibrationMatrixObj.numSamplePerChirp;
numchirps = genCalibrationMatrixObj.nchirp_loops;
numRx = genCalibrationMatrixObj.numRxToEnable;

gBeam_Range_Profile = zeros(numTheta, numFrames, numSamples, numchirps, numRx);
gBeam_RangeDoppler_Profile = zeros(numTheta, numFrames, numSamples, numchirps, numRx);
gBeamGain_power = zeros(1, numTheta);

gBeams_Cube_Angle = zeros(numTheta, numFrames*numchirps, numRx);


% 初始化常用变量
peakValuesBin = zeros(1,genCalibrationMatrixObj.numRxToEnable);
peakValuesTargetDistance = zeros(1, genCalibrationMatrixObj.numRxToEnable);
peakValues = zeros(1, genCalibrationMatrixObj.numRxToEnable);
noiseFloorValues = zeros(1, genCalibrationMatrixObj.numRxToEnable);
peakSNRValues = zeros(1, genCalibrationMatrixObj.numRxToEnable);

% loop through all beam steering angles
for idxTheta = 1:ceil((length(dataFolder)-2)/dataFolder_Skip)
    if (folderIdx > length(dataFolder))
        fprintf(2, "folderIdx is out of the number of dataFolder_calib_data.. %d\n", folderIdx)
        break;
    end

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
    
    for frameId = 1:ceil(numValidFrames/UsedFrame_Skip)
        % use second frame for calibration 
        genCalibrationMatrixObj.frameIdx = frameId; 
        
        % force cal object to single TX phase/frame (numChirpPerLoop <-> numTXPhases)
        %genCalibrationMatrixObj.nchirp_loops = 1;
        genCalibrationMatrixObj.TxToEnable = 1;
        
        % read in .bin files
        % rawData = [samples, chirps, Rx]
        rawData = cascade_Read_rawData(genCalibrationMatrixObj);
        

        % Range FFT  gBeam_Range_Profile(idxTheta, frameId, :, :, :)
        Range_Profile(:, :, :) = fftshift(fft(rawData(:,:,:), genCalibrationMatrixObj.numSamplePerChirp, 1), 1);
        gBeam_Range_Profile(idxTheta,frameId, :, :, :) = Range_Profile(:, :, :);
        

        gBeam_RangeDoppler_Profile(idxTheta, frameId,:,:,:) = fftshift(fft(Range_Profile(:, :, :), ...
                                                                           genCalibrationMatrixObj.nchirp_loops, 2), 2); % * 1/(genCalibrationMatrixObj.nchirp_loops);
        
        % debug plots 
        if(DEBUG_PLOTS)
            plot(axes1, 10*log2(abs(squeeze(Range_Profile( :, :, 1)))));
            title(axes1,"1D Range FFT");
 
            imagesc(axes2, 10*log2(abs(squeeze(gBeam_RangeDoppler_Profile(idxTheta, frameId, :, :, 1)))));
            title(axes2, "After 2D Range-Doppler FFT");           
        end

        for idxRX = 1:genCalibrationMatrixObj.numRxToEnable
            %[TargetBinValue, TargetBinIdx] = max(abs(squeeze(Data_1DFFT(searchBinsSkip:genCalibrationMatrixObj.numSamplePerChirp, genCalibrationMatrixObj.nchirp_loops/2 + 1, 1, idxTheta))));

            Data_2DFFT_Mean_logPower = squeeze(10*log(mean(abs(gBeam_RangeDoppler_Profile(idxTheta, frameId, :, ...
                                                                                            [1:genCalibrationMatrixObj.nchirp_loops/2, ...
                                                                                            genCalibrationMatrixObj.nchirp_loops/2+2:genCalibrationMatrixObj.nchirp_loops], ...
                                                                                            idxRX)), 4)));

            [TargetBinValue, TargetBinIdx] = max(Data_2DFFT_Mean_logPower);

            peakValuesBin(idxRX) = TargetBinIdx;
            offset = 2; % because IF shift need to offset 2
            peakValuesTargetDistance(idxRX) = (TargetBinIdx - (searchBinsSkip - 1) - offset) * genCalibrationMatrixObj.rangeResolution;
          
            % record phase at target peak bin
            peakValues(idxRX) = abs(gBeam_RangeDoppler_Profile(idxTheta, frameId, peakValuesBin(idxRX), genCalibrationMatrixObj.nchirp_loops/2 + 1, idxRX));
            noiseFloorValues(idxRX) = squeeze(mean(abs(gBeam_RangeDoppler_Profile(idxTheta, frameId, peakValuesBin(idxRX), [1:genCalibrationMatrixObj.nchirp_loops/2, genCalibrationMatrixObj.nchirp_loops/2+2:genCalibrationMatrixObj.nchirp_loops],idxRX)), 4)); 
            %peakSNRValues(idxTX, idxRX, idxTheta) = peakValues(idxTX, idxRX, idxTheta) - noiseFloorValues(idxTX, idxRX, idxTheta);
            peakSNRValues(idxRX) = peakValues(idxRX);
            

            % rangeFFT后的增益
            gBeamGain_power(idxTheta) = mean(Data_2DFFT_Mean_logPower(peakValuesBin(idxRX)));
            

            
            % debug plots 
            if(DEBUG_PLOTS)
                plot(axes4, 10*log(abs(squeeze(gBeam_RangeDoppler_Profile(idxTheta, frameId, :, genCalibrationMatrixObj.nchirp_loops/2 + 1, idxRX)))));
                hold(axes4, 'on');
                plot(axes4, Data_2DFFT_Mean_logPower);
                plot(axes4, peakValuesBin(idxRX), Data_2DFFT_Mean_logPower(peakValuesBin(idxRX)), '-o', 'color', 'red');
                hold(axes4, 'off');
                title(axes4, 'Calibration Target IF Spectrum');
                xlabel(axes4, 'Beam Steering Angle (degrees)'); 
                ylabel(axes4, '2D-FFT (0-velocity) Magnitude (dB)');

                plot(axes5, angle(squeeze(gBeam_RangeDoppler_Profile(idxTheta, frameId, :, genCalibrationMatrixObj.nchirp_loops/2 + 1, idxRX))) * 180 / pi );
                hold(axes5, 'on');
                plot(axes5, peakValuesBin(idxRX), angle(squeeze(gBeam_RangeDoppler_Profile(idxTheta, frameId, peakValuesBin(idxRX), genCalibrationMatrixObj.nchirp_loops/2 + 1, idxRX))) * 180 / pi, '-o', 'color', 'red');
                hold(axes5, 'off');                
                title(axes5, 'Calibration Target Phase vs. IF bins');
                xlabel(axes5, '1D-FFT Spectrum (bins)');
                ylabel(axes5, '1D-FFT Phase (degrees)');
                pause(0.01);
            end

        
        end % end for idxRX

        
        
        TargetPeakValueBin = mode(peakValuesBin, 'all');
        TargetpeakValuesDistance = mode(peakValuesTargetDistance, 'all');

        if (frameId == 2)
            disp(["---->>> FrameID=", frameId, "targetPeakValueBin =", TargetPeakValueBin, "rangeResolution =", genCalibrationMatrixObj.rangeResolution, "peakValuesTargetDistance =", TargetpeakValuesDistance]);
        end


    %gBeams_Phase = [gBeams_Phase, Range_Profile(peakValuesBin(1), :, 1)];
    
    end % end for frameId
    
    gPhase_RXs = [];
    for frameId = 1:ceil(numValidFrames/UsedFrame_Skip)
       gPhase_RXs = cat(1, gPhase_RXs, squeeze(gBeam_Range_Profile(idxTheta, frameId, peakValuesBin(1), :, :))); 
    end
    
    % 删除常用局部变量
    clearvars peakValuesBin;
    clearvars peakValuesTargetDistance;
    clearvars peakValues;
    clearvars noiseFloorValues;
    clearvars peakSNRValues;
    clearvars Data_2DFFT_Mean_logPower;
    %
    % ---------------------Processing All Beam Space---------------------
    
    % IQ plan analysis
    if (0)
        figure(6)
        set(gca,'position',[100 100 1500 500]);
        sizeMarker = 10;
        colorMarker = colormap(parula(size(gPhase_RXs, 1)));
        fig6 = tight_subplot(2,8,[.05 .05],[.1 .05],[.05 .05]);
        for i=1:size(gBeams_Cube_Angle, 3)
            axes(fig6(i));
            scatter(real(gPhase_RXs(:, i)), imag(gPhase_RXs(:, i)), sizeMarker, colorMarker);
            colorbar;
            title(sprintf("Rx: %d", i));
        end        
        pause(0.01)
    end
    

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
        % step1 去除 由 number chirp 引起的 相位突变
        angle_gap = 0.0;
        for idxChirp = 1:ceil(numValidFrames/UsedFrame_Skip)-1
            angle_gap = angle_gap + Phase_Angle(idxChirp*genCalibrationMatrixObj.nchirp_loops+1) - Phase_Angle(idxChirp*genCalibrationMatrixObj.nchirp_loops);
            Phase_Angle(idxChirp*genCalibrationMatrixObj.nchirp_loops:(idxChirp+1)*genCalibrationMatrixObj.nchirp_loops-1) = Phase_Angle(idxChirp*genCalibrationMatrixObj.nchirp_loops:(idxChirp+1)*genCalibrationMatrixObj.nchirp_loops-1) - angle_gap;
        end
        

        % step2
        for idxChirp = 1:ceil(numValidFrames/UsedFrame_Skip)-1
            th1 = abs(Phase_Angle(idxChirp*genCalibrationMatrixObj.nchirp_loops-1) -  Phase_Angle(idxChirp*genCalibrationMatrixObj.nchirp_loops));
            th2 = abs(Phase_Angle(idxChirp*genCalibrationMatrixObj.nchirp_loops+1) -  Phase_Angle(idxChirp*genCalibrationMatrixObj.nchirp_loops));
            
            Threshold = abs(Phase_Angle(idxChirp*genCalibrationMatrixObj.nchirp_loops+1)- Phase_Angle(idxChirp*genCalibrationMatrixObj.nchirp_loops-1));
    
            if th1 + th2 > Threshold
                Phase_Angle(idxChirp*genCalibrationMatrixObj.nchirp_loops) = 0.5*(Phase_Angle(idxChirp*genCalibrationMatrixObj.nchirp_loops+1) + Phase_Angle(idxChirp*genCalibrationMatrixObj.nchirp_loops-1));
            end
        end
        
        % step 相位差分
        Beam_Sgnal = diff(Phase_Angle);
        Beam_Sgnal = Beam_Sgnal(1:length(Beam_Sgnal)-1);

        if(DEBUG_PLOTS)
            figure(9);
            plot(unwrap(angle(gPhase_RXs(:, idxRX))),'color', 'blue' );
            hold on;
            plot(Phase_Angle,'color', 'red');
            hold off;
            title("Phase Signal")
            pause(0.01)
        
            % figure(10);
            % plot(Beam_Sgnal);
            % title("gBeams_Phase_Angle Diff")
            % pause(0.01)
        end %if(DEBUG_PLOTS)
            
        gBeams_Cube_Angle(idxTheta, 1:length(Phase_Angle), idxRX) = Phase_Angle;
        
    end % idxRx    
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

    
    %% ------------------ 计算当前 TXBF 角度的 信噪比
    % Range
    



    folderIdx = folderIdx + dataFolder_Skip;
end
% 删除原始数据文件相关的 不用局部变量
clearvars fileIdx_unique;
clearvars fileNameCascade;
clearvars fileNameStruct;
clearvars rawData;


%% 波束增益绘制
figure(11)
subplot(1,2,1)
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
gBeam_RangeDoppler_Profile_zeroDop = squeeze(gBeam_RangeDoppler_Profile(:,2,:,numchirps/2+1,:));
SweepAngles=(1:dataFolder_Skip:(length(dataFolder)-2)) - ceil((length(dataFolder)-2)/2);
range_angle_stich = zeros(size(gBeam_RangeDoppler_Profile_zeroDop, 2), length(SweepAngles));

TI_Cascade_RX_position_azi = [ 11:14 50:53 46:49 0:3  ];
TI_Cascade_RX_ID = [13 14 15 16 1 2 3 4 9 10 11 12 5 6 7 8 ]; %RX channel order on TI 4-chip cascade EVM
D_RX = TI_Cascade_RX_position_azi(TI_Cascade_RX_ID); %RX azimuth antenna coordinate

centerFrequency = Start_Freq_GHz + (genCalibrationMatrixObj.numSamplePerChirp/genCalibrationMatrixObj.Sampling_Rate_sps*Slope_MHzperus)/2;
dlambda = 0.5*centerFrequency/TI_Cascade_Antenna_DesignFreq;

for idxTheta = 1:length(SweepAngles)-1
    angleTX = SweepAngles(idxTheta);
    wx = sind(angleTX);
    a1_az = exp(1j*2*pi*dlambda*(D_RX*wx));
    
    for irange = 1:size(gBeam_RangeDoppler_Profile_zeroDop, 2)
        
        RX_data = squeeze(gBeam_RangeDoppler_Profile_zeroDop(idxTheta, irange, :));
        
        range_angle_stich(irange, idxTheta) = a1_az*(RX_data'*RX_data)*a1_az';
    end
end
sine_theta = sind(SweepAngles);
cos_theta = sqrt(1-sine_theta.^2);
indices_1D_half = 1:size(gBeam_RangeDoppler_Profile_zeroDop, 2)/2;

[R_mat, sine_theta_mat] = meshgrid(indices_1D_half, sine_theta);

[~, cos_theta_mat] = meshgrid(indices_1D_half,cos_theta);

x_axis = R_mat.*cos_theta_mat;
y_axis = R_mat.*sine_theta_mat;

range_angle_stich_half = (range_angle_stich(indices_1D_half + size(gBeam_RangeDoppler_Profile_zeroDop, 2)/2, :).');
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

%%.





%%
TargetBin_Vis = TargetPeakValueBin;
Frame_Vis = 2;

for i = 1:size(gBeam_RangeDoppler_Profile,1)
    figure(13+i)
    imagesc(10*log2(abs(squeeze(gBeam_RangeDoppler_Profile(i, Frame_Vis, :, :, 1)))));
    title("After 2D Range-Doppler FFT");
end


sumDoppler(:,:) = 10*log2(abs(squeeze(sum(gBeam_RangeDoppler_Profile(:, Frame_Vis, :, :, 1),4))));
CM13= colormap(jet(size(sumDoppler,1)));
figure(111)
for i =1:size(sumDoppler,1)
    plot(sumDoppler(i,:), 'Color', CM13(i,:),LineWidth=1);
    hold on;
    leg_str{i} = ['Theta ', num2str(i)];
end
legend(leg_str);

%% 在 Range-Doppler维度找到目标可能存在的峰值

% TODO
% CFAR 检测

cfar_in = gBeam_RangeDoppler_Profile(2, 2, :, :, :);
cfar_model = CFAR_CASO('pfile', pathGenParaFile);
cfar_out = cfar_model.CFAR_RD_Processing(cfar_in);



% 循环对每个Theta角度下峰值检测


% 获取信号能量噪声能量
peakValues(idxRX) = abs(gBeam_RangeDoppler_Profile(:, Frame_Vis, peakValuesBin(idxRX), genCalibrationMatrixObj.nchirp_loops/2 + 1, idxRX));

noiseFloorValues(idxRX) = squeeze(mean(abs(gBeam_RangeDoppler_Profile(idxTheta, frameId, peakValuesBin(idxRX), [1:genCalibrationMatrixObj.nchirp_loops/2, genCalibrationMatrixObj.nchirp_loops/2+2:genCalibrationMatrixObj.nchirp_loops],idxRX)), 4)); 
peakSNRValues(idxTX, idxRX, idxTheta) = peakValues(idxTX, idxRX, idxTheta) - noiseFloorValues(idxTX, idxRX, idxTheta);
peakSNRValues(idxRX) = peakValues(idxRX);
            




fft3D(:,:) = fftshift(fft(gBeam_Range_Profile(10, 2, :, numchirps/2+1, :), 180, 5), 5);
figure(122)
imagesc(10*log10(abs(squeeze(fft3D))));
title("After 3D Range-Angle FFT");
pause(0.01);


%% 绘制所有beam 的 Range-Angle 图像, 利用 RA图 确定目标的位置角度


% gBeam_RangeAngle_Profile(idxTheta, frameId,:,:,:)

% doppler 维度求均值
% temp = mean(gBeam_RangeAngle_Profile, 4);




%%  在 Beam 所扫描的平面上寻找最优区域，计算 PSINR

% step1. 计算beam噪声平面,在RX1
nFFT = 1024;
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








toc
timeProcessing = toc;


