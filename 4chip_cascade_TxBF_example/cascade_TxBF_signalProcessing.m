%  Copyright (C) 2018 Texas Instruments Incorporated - http://www.ti.com/
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

% read the adc data with advanced frame configuration with TX beamforming,
% stich the multiple scans together with RX beamforming

%read the adc data and postprocess it
close all;clc
clearvars;
addpath(genpath('.\PostProc\'))
addpath(genpath('.\cascade_json_parser\'))

%load the adc data params
%BF_data_dest0 = 'E:\ti\mmwave_studio_03_00_00_14\mmWaveStudio\PostProc\TXBF_BeamAngle31\';
%% todo
BF_data_dest0 = 'E:\ti\mmwave_studio_03_00_00_14\mmWaveStudio\PostProc\TxBF_USRR_BeamAngle_test2\';

%% define what angles to steer in TX beamforming mode
paramsConfig.anglesToSteer=[1]; 

% Only Frame based supported in the Post-processing script
paramsConfig.Chirp_Frame_BF = 0; % 1 - chirp based beam steering, 0 - frame based beam steering

paramsConfig.NumAnglesToSweep = length(paramsConfig.anglesToSteer);

%parse Json file
paramsConfig = parameter_gen_from_Jason(BF_data_dest0, paramsConfig);

%load MIMO calbibration params
calibrationFilePath = 'PostProc\calibrateResults_high.mat';
load(calibrationFilePath);
paramsConfig.Slope_MHzperus_calib = params.Slope_MHzperus;
BF_MIMO_ref = calibResult.RxMismatch;

%re-arrange calibration data
RangeMat = calibResult.RangeMat;
RangeMat = RangeMat(:, paramsConfig.TI_Cascade_RX_ID);
PeakValMat = calibResult.PeakValMat;
PeakValMat = PeakValMat(:, paramsConfig.TI_Cascade_RX_ID);

%pass parameters from paramsConfig
numFrames = paramsConfig.Num_Frames;
Rx_Ant_Arr = paramsConfig.TI_Cascade_RX_ID;
anglesToSteer = paramsConfig.anglesToSteer;

currentFolder = pwd;
cd(BF_data_dest0);
listing = dir('*_data.bin');
cd(currentFolder);
cnt = 1;

   % Get Unique File Idxs in the "dataFolder_test"   
   [fileIdx_unique] = getUniqueFileIdx(BF_data_dest0)
  
  
for i_file = 1:(length(fileIdx_unique))
    
     % Get File Names for the Master, Slave1, Slave2, Slave3   
   [fileNameStruct]= getBinFileNames_withIdx(BF_data_dest0, fileIdx_unique{i_file});

    range_angle_stich = [];
    RX_TXBF_calVec = ones(1,length(paramsConfig.Rx_Elements_To_Capture));
    
    % Get Valid Number of Frames 
    [numValidFrames, dataFileSize] = getValidNumFrames(fullfile(BF_data_dest0, fileNameStruct.masterIdxFile));
    
    % beam 处理数据需要的变量
    gBeams_Phase = []; %zeros(length(paramsConfig.anglesToSteer),1)
    gDownsample_value = 2;
    
    doppler_fs = 1e6 / (paramsConfig.SF1SubFramePeriodicity / 200) * paramsConfig.SF1NumLoops; 


    for frameId = 2:numValidFrames
        paramsConfig.frameId = frameId;
        [radar_data_TXBF]= fcn_read_AdvFrmConfig_BF_Json(fileNameStruct,paramsConfig);
      
        % TI_Cascade_RX_ID: RX channel ID order to be read in sequence
        % radar_data_TXBF[nSamples, nChirps, nRx, nBeams]
        radar_data_TXBF = radar_data_TXBF(:,:,paramsConfig.TI_Cascade_RX_ID,:);
        % -----------------process TXBF data
        if frameId < 1
            figure(1)
            startRangeInd = 5;
            lastRangeIndToThrow = 20;
            
            [ range_angle_stich]= Plot_advFraConfig_TXBF_rangeAzimuth_stich( ...
                                        radar_data_TXBF,...
                                        Rx_Ant_Arr, ...
                                        paramsConfig, ...
                                        anglesToSteer, ...
                                        BF_MIMO_ref(Rx_Ant_Arr), ...
                                        startRangeInd, ...
                                        lastRangeIndToThrow);
            
            pause(0.1)

        end
        
        % radar_data_TXBF
        rangeFFT = fft(radar_data_TXBF, paramsConfig.Samples_per_Chirp, 1);
        DopplerFFT = fft(rangeFFT, paramsConfig.nchirp_loops, 2);            

        % estimate the range bin, plot all beam angle range fft         
        if frameId < 100
            figure(12);
            for i = 1:length(paramsConfig.anglesToSteer)
                plot(abs(rangeFFT(:, 1, 1, i)), LineWidth=1);
                hold on;
            end
            title("Range FFT")
            pause(0.05);
        end 

        % 获取 不同 beam 的信号数据
        Rangebin = 9;  % 人所在的位置
        for i = 1:length(paramsConfig.anglesToSteer)
            beam_phase = rangeFFT(Rangebin, :, 1, i);
            gBeams_Phase = [gBeams_Phase, beam_phase];
        end                
    end
    cnt = cnt + 1;  % 用于记录不同的次数
    
    %%
    % 所有 beam 的 数据
    gBeams_All_Span = gBeams_Phase

    gBeams_Phase = reshape(gBeams_Phase, ...
                       length(paramsConfig.anglesToSteer), ...
                       int32(length(gBeams_Phase)/length(paramsConfig.anglesToSteer)));

    
    % 1. 相位变化一个复数， a+j.b  angle-> (b/a) ==> unwrap() == > diff 
    % 可视化数据，将每个 beam 的数据可视化
    figure(100);
    for i = 1:length(paramsConfig.anglesToSteer)
        gBeams_motion(i,:) = diff(unwrap(angle(gBeams_Phase(i,:))));

        x = ((i-1)*length(gBeams_motion(i,:))+1:i*length(gBeams_motion(i,:)));
        plot(x, gBeams_motion(i,:));
        hold on;
        pause(0.01);
    end
    
    % 低通滤波？
    % 2. filter 带通滤波 呼吸：0.3 - 0.8； 心跳频率 0.8 - 2Hz
    h = ellip_band_filter;
    figure(101);
    for i = 1:length(paramsConfig.anglesToSteer)
        gBeams_motion_filter(i,:) = filter(h, gBeams_motion(i,:));
        x = ((i-1)*length(gBeams_motion(i,:))+1:i*length(gBeams_motion(i,:)));
        plot(x, gBeams_motion_filter(i,:));
        hold on;
        pause(0.01);
    end
    
    % 峰值检测
    % pks = findpeaks(gBeams_motion_filter(0,:), MinPeakDistance=doppler_fs/2);

end
%%
