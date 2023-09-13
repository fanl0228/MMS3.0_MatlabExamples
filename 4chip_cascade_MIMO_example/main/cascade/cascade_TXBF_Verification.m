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
clearvars
close all

pro_path = getenv('CASCADE_SIGNAL_PROCESSING_CHAIN_MIMO');
input_path = strcat(pro_path,'\main\cascade\input\');
dataPlatform = 'TDA2'

DEBUG_PLOTS = 1;                % optionally display debug plots while calibration data is being processed 

numDevices = 4;                 % number of AWRx devices being processed, set to 4, for the full MMWCAS-RF-EVM device array
numTX = 3;                      % number of AWRx TX channels being processed, set to 3, for the full AWRx device channels
numRX = 16;                     % number of MMWCAS-RF-EVM RX channels being processed, set to 16, for the full MMWCAS-RF-EVM RX channels
numPhaseShifterOffsets = 64;    % number of phase-shifter offset increments being processed, set to 64 for full phase-shifter range (number of datasets)
numChirpsLoopsPerFrame = 12;    % number of chirp-loops per frame. Only a single TX active per chirp-loop. 
numChirpsPerLoop = 128;         % number of chirp loops
numSamplesPerChirp = 256;       % number of samples per chirp
searchBinsSkip = 40;            % number of bins to skip when looking for peak from corner reflector 

% select reference TX/Device channel for computing offsets - determined by 
% TX antenna array geometry and phase-shifter offset utilization (all RX
% channels are processed)
refTX = 10;
refPhaseOffset = 1;

targetRange = 2.3;  % estimated corner-reflector target range (in meters) 
dataFolder_calib_data_path = 'C:\Radar_Data\20200721_phase_shifter_cal_testing\20200723_testdata_SN5733600017_outdoor_TXBF_test1\nocal\'; % folder holding all of the calibration datasets
%dataFolder_calib_data_path = 'C:\Radar_Data\20200721_phase_shifter_cal_testing\20200723_testdata_SN5733600017_outdoor_TXBF_test1\cal\'; % folder holding all of the calibration datasets

fig1 = figure(1); % FFT magnitude (dB)
axes1 = axes;

fig2 = figure(2); % FFT phase (deg)
axes2 = axes;

fig3 = figure(3); % accumulated target bin values
axes3 = axes;

fig4 = figure(4); % accumulated target distance values
axes4 = axes;

fig5 = figure(5); % accumulated target phase angle values
axes5 = axes;


%parameter file name for the test
pathGenParaFile = [input_path,'generateClibrationMatrix_param.m'];

%important to clear the same.m file, since Matlab does not clear cache
%automatically
clear(pathGenParaFile);



% loop through each folder in the calibration directory create 
% matrix with dimensions [devices, number of TX, number of phase shift offsets] 
% with the phase values for the detected range-bin peak
dataFolder_calib_data_info = dir(dataFolder_calib_data_path);

iterationMax = numPhaseShifterOffsets;
iteration = 0;

% find the first folder of cal data
for folderIdx = 1:length(dataFolder_calib_data_info)
    if(dataFolder_calib_data_info(folderIdx).isdir && dataFolder_calib_data_info(folderIdx).name ~= "." && dataFolder_calib_data_info(folderIdx).name ~= "..") 
        folderIdxStart = folderIdx;
        break;
    end

end

folderIdx = folderIdxStart; % start at first index that is a folder

for idxTheta = 1:31 % loop through all beam steering angles

    if(dataFolder_calib_data_info(folderIdx).isdir) 
        disp(['dataFolder_calib_data_info(folderIdx).name = ', dataFolder_calib_data_info(folderIdx).name]);

        % create folder structure for this iteration        
        fileNameCascade.dataFolderName = strcat([dataFolder_calib_data_info(folderIdx).folder, '\', dataFolder_calib_data_info(folderIdx).name, '\']);
        [fileIdx_unique] = getUniqueFileIdx(fileNameCascade.dataFolderName);
        [fileNameStruct] = getBinFileNames_withIdx(fileNameCascade.dataFolderName, fileIdx_unique{1});                  

        %generate parameter file for the dataset
        parameter_file_gen_TXPS_VerificationCalib_json(fileNameCascade.dataFolderName, pathGenParaFile, dataPlatform);

        %generate calibration matrix object for the dataset
        genCalibrationMatrixObj = genCalibrationMatrixCascade('pfile', pathGenParaFile,...
        'calibrateFileName',fileNameCascade.dataFolderName, 'targetRange', targetRange);
        genCalibrationMatrixObj.binDataFile = fileNameStruct;% dataFolder_calib_data;%[dataFolder_calib_data listing.name];
    end
            
    % use second frame for calibration 
    genCalibrationMatrixObj.frameIdx = 2;             
    
    % force cal object to single TX phase/frame (numChirpPerLoop <-> numTXPhases)
    %genCalibrationMatrixObj.nchirp_loops = 1;
    genCalibrationMatrixObj.TxToEnable = 1;
    
    % read in .bin files
    %calData(deviceIdx, TXIdx, PSIdx, :, :, :) = cascade_Read_TX_Cal_Data(genCalibrationMatrixObj);
    calData = cascade_Read_TX_Cal_Data(genCalibrationMatrixObj);
          

    for idxTX = 1:1 % loop through each TX phase   
        for idxRX = 1:1 % loop through each RX        

            calData_1DFFT(:, :, idxTX, idxRX, idxTheta) = fftshift(fft(calData(:, :, idxRX, idxTX), genCalibrationMatrixObj.numSamplePerChirp, 1));
            calData_2DFFT(:, :, idxTX, idxRX, idxTheta) = 1/(numChirpsPerLoop) * fftshift(fft(calData_1DFFT(:, :, idxTX, idxRX, idxTheta), numChirpsPerLoop, 2));
            
          
            % find target peak bin in 2D-FFT 0-velocity bin (skip close
            % bins to avoid DC leakage or bumper reflections
            [TargetBinValue, TargetBinIdx] = max(abs(squeeze(calData_2DFFT(searchBinsSkip:numSamplesPerChirp, numChirpsPerLoop/2 + 1, idxTX, idxRX, idxTheta))));
            peakValuesBin(idxTX, idxRX, idxTheta) = TargetBinIdx + searchBinsSkip - 1;
            %peakValuesBin(idxTX, idxRX, idxTheta) = TargetBinIdx;
            peakValuesTargetDistance(idxTX, idxRX, idxTheta) = TargetBinIdx * genCalibrationMatrixObj.rangeResolution;

            % record phase at target peak bin
            peakValues(idxTX, idxRX, idxTheta) = abs(calData_2DFFT(peakValuesBin(idxTX, idxRX, idxTheta), numChirpsPerLoop/2 + 1, idxTX, idxRX, idxTheta));
            noiseFloorValues(idxTX, idxRX, idxTheta) = mean(abs(squeeze(calData_2DFFT(peakValuesBin(idxTX, idxRX, idxTheta), [1:numChirpsPerLoop/2, numChirpsPerLoop/2+2:numChirpsPerLoop], idxTX, idxRX, idxTheta))), 2); 
            %peakSNRValues(idxTX, idxRX, idxTheta) = peakValues(idxTX, idxRX, idxTheta) - noiseFloorValues(idxTX, idxRX, idxTheta);
            peakSNRValues(idxTX, idxRX, idxTheta) = peakValues(idxTX, idxRX, idxTheta);
                               
            % debug plots 
            if(DEBUG_PLOTS)


                plot(axes1, 10*log(abs(calData_2DFFT(:, numChirpsPerLoop/2 + 1, idxTX, idxRX, idxTheta))));
                hold(axes1, 'on');
                plot(axes1, 10*log(mean(abs(squeeze(calData_2DFFT(:, [1:numChirpsPerLoop/2, numChirpsPerLoop/2+2:numChirpsPerLoop], idxTX, idxRX, idxTheta))), 2)));
                plot(axes1, peakValuesBin(idxTX, idxRX, idxTheta), 10*log(peakValues(idxTX, idxRX, idxTheta)), '-o', 'color', 'red');
                hold(axes1, 'off');
                title(axes1, 'Calibration Target IF Spectrum');
                xlabel(axes1, 'Beam Steering Angle (degrees)'); 
                ylabel(axes1, '2D-FFT (0-velocity) Magnitude (dB)');

                plot(axes3, squeeze(peakValuesBin(idxTX, idxRX, :)));                
                title(axes3, 'Calibration Target Detected Index');
                xlabel(axes3, 'Beam Steering Angle (degrees)');
                ylabel(axes3, 'Calibration Target Sampled IF Index');


                plot(axes4, squeeze(peakValuesTargetDistance(idxTX, idxRX, :)));
                title(axes4, 'Calibration Target Distance');
                xlabel(axes4, 'Beam Steering Angle (degrees)');
                ylabel(axes4, 'Target Distance (meters)');


                plot(axes5, squeeze(10*log(peakValues(idxTX, idxRX, :))));
                title(axes5, 'Calibration Target Phase vs. IF bins');
                xlabel(axes5, 'Beam Steering Angle (degrees)');
                ylabel(axes5, 'Peak 2D-FFT (0-velocity) Magnitude (dB)');            
                
                pause(0.01);


            end
            
            
        end   
        
        
        
    end

    folderIdx = folderIdx + 1;

end


toc
timeProcessing = toc;


