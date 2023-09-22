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

% cascade_TX_Phase_Calibration.m
%  
% Top level main test chain to perform antenna calibration. 
% genCalibrationMatrixCascade module is first initialized before
% actually used in the chain. The output is saved to a .mat file to be used
% later in data processing

tic
% clearvars
% close all
setenv('CASCADE_SIGNAL_PROCESSING_CHAIN_MIMO', 'D:\ti\mmwave_studio_03_00_00_14\mmWaveStudio\MMS3.0_MatlabExamples\4chip_cascade_MIMO_example')

pro_path = getenv('CASCADE_SIGNAL_PROCESSING_CHAIN_MIMO');
input_path = strcat(pro_path,'\main\cascade\input\');
dataPlatform = 'TDA2'

DEBUG_PLOTS = 0;                % optionally display debug plots while calibration data is being processed 

searchBinsSkip = 129;            % number of bins to skip when looking for peak from corner reflector 

targetRange = 2.0;  % estimated corner-reflector target range (in meters) 

BeamAngle_Data_Path = 'I:\20230915_TXBF_AngleSweep_Vib_Range2.1m_test3\Test_Vib_TXBF_BeamAngle01\'; 
% folder holding all of the calibration datasets

if(DEBUG_PLOTS) 
    fig1 = figure(1); % FFT magnitude (dB)
    axes1 = axes;
    
    fig2 = figure(2); % FFT phase (deg)
    axes2 = axes;
    
end
datacolor = colormap(jet(20));



%parameter file name for the test
pathGenParaFile = [input_path,'generateClibrationMatrix_param.m'];

%important to clear the same.m file, since Matlab does not clear cache
%automatically
clear(pathGenParaFile);



%iterationMax = numPhaseShifterOffsets;
iteration = 0;
gBeamGain = [];

[fileIdx_unique] = getUniqueFileIdx(BeamAngle_Data_Path);

[fileNameStruct] = getBinFileNames_withIdx(BeamAngle_Data_Path, fileIdx_unique{1});

[numValidFrames, dataFileSize] = getValidNumFrames(fullfile(BeamAngle_Data_Path, fileNameStruct.masterIdxFile));

%generate parameter file for the dataset
parameter_AngleSweep_json(BeamAngle_Data_Path, pathGenParaFile, dataPlatform);

%generate calibration matrix object for the dataset
genCalibrationMatrixObj = genCalibrationMatrixCascade('pfile', pathGenParaFile,...
'calibrateFileName',BeamAngle_Data_Path, 'targetRange', targetRange);
genCalibrationMatrixObj.binDataFile = fileNameStruct;% dataFolder_calib_data;%[dataFolder_calib_data listing.name];


gBeams_Phase = []

for frameId = 1:numValidFrames
    % use second frame for calibration 
    genCalibrationMatrixObj.frameIdx = frameId; 
    
    % force cal object to single TX phase/frame (numChirpPerLoop <-> numTXPhases)
    %genCalibrationMatrixObj.nchirp_loops = 1;
    genCalibrationMatrixObj.TxToEnable = 1;
    
    % read in .bin files
    %calData(deviceIdx, TXIdx, PSIdx, :, :, :) = cascade_Read_TX_Cal_Data(genCalibrationMatrixObj);
    rawData = cascade_Read_rawData(genCalibrationMatrixObj);

    % Range FFT
    Data_1DFFT(:,:,:) = fftshift(fft(rawData(:,:,:), genCalibrationMatrixObj.numSamplePerChirp, 1), 1);
    Data_2DFFT(:,:,:) = fftshift(fft(Data_1DFFT(:,:,:), genCalibrationMatrixObj.nchirp_loops, 2),2); % * 1/(genCalibrationMatrixObj.nchirp_loops);
    Data_3DFFT(:,:,:) = fftshift(fft(Data_1DFFT(:,:,:), 180, 3), 3);
    
    % debug plots 
    if(DEBUG_PLOTS)
        figure(100)
        for i = 1:12
            plot(10*log2(abs(Data_1DFFT(:, 1, i))));
            hold on;
        end
        hold off;
        title("1D Range FFT");
        pause(0.05);
                        
        figure(101)
        imagesc(10*log2(abs(Data_2DFFT(:, :, 1))));
        title("After 2D Range-Doppler FFT");
        pause(0.05);

        figure(102)
        Data_3DFFT_vis = permute(Data_3DFFT,[1,3,2,4]);
        imagesc(10*log2(abs(Data_3DFFT_vis(:, :, 1))));
        title("After 3D Range-Angle FFT");
        pause(0.05);

    end
    
    for idxRX = 1:16
        %[TargetBinValue, TargetBinIdx] = max(abs(squeeze(Data_1DFFT(searchBinsSkip:genCalibrationMatrixObj.numSamplePerChirp, genCalibrationMatrixObj.nchirp_loops/2 + 1, 1))));
        
        Data_2DFFT_Mean_logPower = 10*log(mean(abs(squeeze(Data_2DFFT(:, [1:genCalibrationMatrixObj.nchirp_loops/2, genCalibrationMatrixObj.nchirp_loops/2+2:genCalibrationMatrixObj.nchirp_loops], idxRX))), 2));
        
        [TargetBinValue, TargetBinIdx] = max(Data_2DFFT_Mean_logPower);

        peakValuesBin(idxRX) = TargetBinIdx; %+ searchBinsSkip - 1;
        
        peakValuesTargetDistance(idxRX) = (TargetBinIdx - (searchBinsSkip - 1)) * genCalibrationMatrixObj.rangeResolution;

        % record phase at target peak bin
        peakValues(idxRX) = abs(Data_2DFFT(peakValuesBin(idxRX), genCalibrationMatrixObj.nchirp_loops/2 + 1, idxRX));
        noiseFloorValues(idxRX) = mean(abs(squeeze(Data_2DFFT(peakValuesBin(idxRX), [1:genCalibrationMatrixObj.nchirp_loops/2, genCalibrationMatrixObj.nchirp_loops/2+2:genCalibrationMatrixObj.nchirp_loops],idxRX))), 2); 
        %peakSNRValues(idxTX, idxRX) = peakValues(idxTX, idxRX) - noiseFloorValues(idxTX, idxRX);
        peakSNRValues(idxRX) = peakValues(idxRX);
        
        % rangeFFT后的增益
        gBeamGain(frameId, idxRX) = Data_2DFFT_Mean_logPower(peakValuesBin(idxRX));

        % debug plots 
        if(DEBUG_PLOTS)
            plot(axes1, 10*log(abs(Data_2DFFT(:, genCalibrationMatrixObj.nchirp_loops/2 + 1, idxRX))));
            hold(axes1, 'on');
            plot(axes1, Data_2DFFT_Mean_logPower);
            plot(axes1, peakValuesBin(idxRX), Data_2DFFT_Mean_logPower(peakValuesBin(idxRX)), '-o', 'color', 'red');
            hold(axes1, 'off');
            title(axes1, 'Calibration Target IF Spectrum');
            xlabel(axes1, 'Beam Steering Angle (degrees)'); 
            ylabel(axes1, '2D-FFT (0-velocity) Magnitude (dB)');

            plot(axes2, angle(squeeze(Data_2DFFT(:, genCalibrationMatrixObj.nchirp_loops/2 + 1,idxRX))) * 180 / pi );
            hold(axes2, 'on');
            plot(axes2, peakValuesBin(idxRX), angle(squeeze(Data_2DFFT(idxRX))) * 180 / pi, '-o', 'color', 'red');
            hold(axes2, 'off');                
            title(axes2, 'Calibration Target Phase vs. IF bins');
            xlabel(axes2, '1D-FFT Spectrum (bins)');
            ylabel(axes2, '1D-FFT Phase (degrees)');

         

            pause(0.01);
        end

    end % end for idxRX
    
    
    gBeams_Phase = [gBeams_Phase, Data_1DFFT(peakValuesBin(idxRX), :, 1)];
end % end for frameId

if (1)
    figure(103);
    scatter(real(gBeams_Phase), imag(gBeams_Phase));
    pause(0.01)
    
    
    gBeams_Phase_Angle = unwrap(angle(gBeams_Phase));
    figure(104);
    plot(gBeams_Phase_Angle);
    title("gBeams_Phase_Angle = unwrap(angle(gBeams_Phase))")
    pause(0.01)

    if (0)
        gBeam_Sgnal1 = diff(gBeams_Phase_Angle);
        gBeam_Sgnal1 = gBeam_Sgnal1(:,1:length(gBeam_Sgnal1)-1);
        figure(400);
        plot(gBeam_Sgnal1);
        pause(0.01)
    end


    % 去除 由 number chirp 引起的 相位突变
    angle_gap = 0.0;
    for idxChirp = 1:numValidFrames-1
        angle_gap = angle_gap + gBeams_Phase_Angle(:, idxChirp*genCalibrationMatrixObj.nchirp_loops+1) - gBeams_Phase_Angle(:, idxChirp*genCalibrationMatrixObj.nchirp_loops);
        gBeams_Phase_Angle(:, idxChirp*genCalibrationMatrixObj.nchirp_loops:(idxChirp+1)*genCalibrationMatrixObj.nchirp_loops-1) = gBeams_Phase_Angle(:, idxChirp*genCalibrationMatrixObj.nchirp_loops:(idxChirp+1)*genCalibrationMatrixObj.nchirp_loops-1) - angle_gap;
    end
    hold on;
    plot(gBeams_Phase_Angle,'color', 'k');


    for idxChirp = 1:numValidFrames-1
        th1 = abs(gBeams_Phase_Angle(:, idxChirp*genCalibrationMatrixObj.nchirp_loops-1) -  gBeams_Phase_Angle(:, idxChirp*genCalibrationMatrixObj.nchirp_loops));
        th2 = abs(gBeams_Phase_Angle(:, idxChirp*genCalibrationMatrixObj.nchirp_loops+1) -  gBeams_Phase_Angle(:, idxChirp*genCalibrationMatrixObj.nchirp_loops));
        
        Threshold = abs(gBeams_Phase_Angle(:, idxChirp*genCalibrationMatrixObj.nchirp_loops+1)- gBeams_Phase_Angle(:, idxChirp*genCalibrationMatrixObj.nchirp_loops-1));

        if th1 + th2 > Threshold
            gBeams_Phase_Angle(:, idxChirp*genCalibrationMatrixObj.nchirp_loops) = 0.5*(gBeams_Phase_Angle(:, idxChirp*genCalibrationMatrixObj.nchirp_loops+1) + gBeams_Phase_Angle(:, idxChirp*genCalibrationMatrixObj.nchirp_loops-1));
        end
    end
    hold on;
    plot(gBeams_Phase_Angle,'color', 'r');
    hold off;
    pause(0.01)


    figure(105);
    plot(gBeams_Phase_Angle);
    title("gBeams_Phase_Angle after phase smooth")
    pause(0.01)
    
    gBeam_Signal = diff(gBeams_Phase_Angle);
    gBeam_Signal = gBeam_Signal(:,1:length(gBeam_Signal)-1);
    
    figure(107);
    plot(gBeam_Signal);
    title("gBeams_Phase_Angle diff")
    pause(0.01)

    
    % 计算振动频率
    signal_fft = 10*log(abs(fftshift(fft(gBeam_Signal, 1024))));
    figure(108)
    plot(signal_fft)
    pause(0.01)


end




toc
timeProcessing = toc;


