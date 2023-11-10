%
%
%
%

function [Intensity_estSINR, Phase_estSINR, gframe_obj] ...
            = Each_Steering_Calculate_pSINR(dataFolder_Path, ...
                                            isStepAngle, ...
                                            TxBF_Angle, ...
                                            RxBF_Angle, ...
                                            UsedFrames, ...
                                            NumRx, ...
                                            num_chirps, ...
                                            Signal_FS,...
                                            LOG_ON, ...
                                            PLOT_ON, ...
                                            LogFileId)

if (isStepAngle == 1)
    dataFolderName = dataFolder_Path;
else
    dataFolder = dir(dataFolder_Path);
    floder_offset = 2;
    folderIdx = TxBF_Angle + ceil((length(dataFolder) - 2) / 2) + floder_offset;
    dataFolderName = strcat([dataFolder(folderIdx).folder, '\', dataFolder(folderIdx).name, '\']);
end

gframe_obj = Each_Steering_processing(dataFolderName, ...
                                        TxBF_Angle, ...
                                        RxBF_Angle, ...
                                        UsedFrames, ...
                                        LOG_ON, ...
                                        PLOT_ON);


%% ---------------------- connect all chirps, and extract phase processing 
gPhase_RXs = [];
Intensity_SNR = [];
targetbins = [];
for frameId = 1:length(gframe_obj)
    targetbins(frameId) = gframe_obj{frameId}.rangeBins(1);
end
% 找 所有frame 中 target的众数
targetbin = mode(targetbins);

for frameId = 1:length(gframe_obj)

    % targetbin = gframe_obj{frameId}.rangeBins(1);

    % Calculate first part Sensing Signal Intensity SNR
    target_idx = find(gframe_obj{frameId}.rangeBins == targetbin);
    if (isempty(target_idx))
        if LOG_ON
            disp("skip ======>>>>> FrameId = ", num2str(frameId));
        end
        continue;
    end

    locs = setdiff(1:length(gframe_obj{frameId}.RA_Spectrum), targetbin+1);
    Intensity_RA_SNR(frameId) = 10*log10((abs(gframe_obj{frameId}.RA_Spectrum(targetbin+1)).^2)./abs(mean(gframe_obj{frameId}.RA_Spectrum(locs))));

    Intensity_SNR(frameId) = sum([gframe_obj{frameId}.estSNR(target_idx)]) / length(target_idx);

    Intensity_SNR(frameId) =  Intensity_SNR(frameId) + Intensity_RA_SNR(frameId);

    if LOG_ON
        disp({strcat("======>>>>> FrameId = ", num2str(frameId)), ...
              strcat("Target Range bin = ", num2str(targetbin)), ...
              strcat("Objects Id= ", num2str(target_idx)), ...
              strcat("Sub_intensity_SNR(dB) = ", num2str(Intensity_SNR(frameId)))});

        fprintf(LogFileId, '%s,\t %s,\t %s,\t %s \n', ...
                        strcat("======>>>>> FrameId = ", num2str(frameId)), ...
                        strcat("Target Range bin = ", num2str(targetbin)), ...
                        strcat("Objects Id= ", num2str(target_idx)), ...
                        strcat("Sub_intensity_SNR(dB) = ", num2str(Intensity_SNR(frameId))));

    end

    % Concatenate the phase data of each frame
    if (targetbin == 0)
        continue;
    else
        gRange_Profile = gframe_obj{frameId}.Range_Profile;
        gPhase_RXs = cat(1, gPhase_RXs, squeeze(gRange_Profile(targetbin, :, :) ) );
    end
end

% --------------- Calculate the intensity signal-to-noise ratio
if isempty(Intensity_SNR)
    Intensity_estSINR = -1;
    Phase_estSINR = -1;
    return;
else
    Intensity_estSINR = sum(Intensity_SNR) / length(Intensity_SNR);
end
if LOG_ON
    disp(strcat("======>>>>> Intensity_SNR(dB) = ", num2str(Intensity_estSINR)));

    fprintf(LogFileId, '%s \n', strcat("======>>>>> Intensity_SNR(dB) = ", num2str(Intensity_estSINR)));
end

% IQ plan analysis the phase data
if 0
    figure(130)
    set(gcf,'units','normalized','outerposition',[0.1 0.1 0.8 0.8])
    sizeMarker = 10;
    colorMarker = colormap(parula(size(gPhase_RXs, 1)));
    fig6 = tight_subplot(2,4,[.05 .05],[.1 .05],[.05 .05]);
    for i=1:size(gPhase_RXs, 2)/2
        axes(fig6(i));
        scatter(real(gPhase_RXs(:, i)), imag(gPhase_RXs(:, i)), sizeMarker, colorMarker);
        colorbar;
        title(sprintf("Rx: %d", i));
    end        
    pause(0.01)
end

%%
% The phase signal is calibrated to facilitate the extraction of fine-grained motion signals
gCalib_Angle = frame_phase_calibration(gPhase_RXs, length(gPhase_RXs)/num_chirps, NumRx, num_chirps, PLOT_ON);

% Calculate second part Sensing Signal Intensity SNR
Phase_estSINR = Calculate_Motion_SINR(gCalib_Angle, Signal_FS, length(gCalib_Angle)/num_chirps, num_chirps, NumRx);

est_PSINR = Intensity_estSINR + Phase_estSINR;

if LOG_ON
    disp(["Intensity_estSINR(dB): ", Intensity_estSINR, ...
            " Phase_estSINR(dB): ", Phase_estSINR, ...
            " est_PSINR(dB): ", est_PSINR])

    fprintf(LogFileId, '%s, \t %s, \t %s \n', ...
                    "Intensity_estSINR(dB): ", Intensity_estSINR, ...
                    " Phase_estSINR(dB): ", Phase_estSINR, ...
                    " est_PSINR(dB): ", est_PSINR);
    fprintf(LogFileId, '\n\n');
end

%% 可视化 差分信号
if PLOT_ON
      % All RX antennas are averaged and the calibrated data signal is plotted.
    gAngle_tmp = mean(unwrap(angle(gPhase_RXs(:, :))), 2); 
    gAngle_tmp_diff = diff(gAngle_tmp);
    
    % ---------------------Plot figure
    fig8=figure(130);
    %set(gcf,'units','normalized','outerposition', [0.7 0.6 0.3 0.4]);
    subplot(121)
    plot(gAngle_tmp,'color', 'blue');
    title({"Phase", strcat("Target Bin",num2str(targetbin))});
    subplot(122)
    plot(gAngle_tmp_diff,'color', 'blue');
    title({"No phase CAL Motion Signal", ...
            strcat("Intensity\_estSINR (dB): ", num2str(Intensity_estSINR)),...
            strcat("Phase\_estSINR (dB): ", num2str(Phase_estSINR))});
    ylim([-1, 1]);
    pause(0.01)


    % All RX antennas are averaged and the calibrated data signal is plotted.
    gAngle = mean(gCalib_Angle(:, :), 2);
    gAngle_diff = diff(gAngle);
    
    % ---------------------Plot figure
    fig9=figure(131);
    %set(gcf,'units','normalized','outerposition', [0.7 0.6 0.3 0.4]);
    subplot(121)
    plot(gAngle,'color', 'blue');
    title({"Phase", strcat("Target Bin",num2str(targetbin))});
    subplot(122)
    plot(gAngle_diff,'color', 'blue');
    title({"Motion Signal", ...
            strcat("Intensity\_estSINR (dB): ", num2str(Intensity_estSINR)),...
            strcat("Phase\_estSINR (dB): ", num2str(Phase_estSINR))});
    ylim([-1, 1]);
    pause(0.01)

    %---------------------- Save figure
    temp = split(dataFolderName, '\');
    angleName = temp(end-1);
    temp(end)=[];
    temp(end)=[];
    tempStr = temp(1);
    for iTemp = 2:length(temp)
        tempStr = strcat(tempStr,'\', temp(iTemp));
    end
    tempStr = cell2mat(tempStr);
    png_floder = strcat(tempStr, '_png\');
    if ~exist(png_floder,'dir')
        mkdir(png_floder)
    end
    png_file = strcat([png_floder, cell2mat(angleName), '.png']);
    saveas(fig9, png_file, 'png');
    fig_file = strcat([png_floder, cell2mat(angleName), 'fig']);
    saveas(fig9, fig_file, 'fig');

end

end

