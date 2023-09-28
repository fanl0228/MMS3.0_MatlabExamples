%
%
%
%

function [est_PSINR, Txbeam_angle] = Each_Steering_Calculate_pSINR(dataFolder_Path, Txbeam_angle, UsedFrames, NumRx, num_chirps, LOG_ON, PLOT_ON)
 
dataFolder = dir(dataFolder_Path);
floder_offset = 2;
folderIdx = Txbeam_angle + ceil((length(dataFolder) - 2) / 2) + floder_offset;
dataFolderName = strcat([dataFolder(folderIdx).folder, '\', dataFolder(folderIdx).name, '\']);


[gframe_obj, gRange_Profile, gRangeDoppler_Profile] = Each_Steering_processing(dataFolderName, Txbeam_angle, UsedFrames, LOG_ON, PLOT_ON);


%% connect all chirps, and extract phase processing 
gPhase_RXs = [];
for frameId = 1:length(gframe_obj)

    targetbin = gframe_obj{frameId}.rangeBins(1);
    
    % Calculate first part Sensing Signal Intensity SNR
    target_idx = find(gframe_obj{frameId}.rangeBins == targetbin);
    Intensity_SNR(frameId) = sum([gframe_obj{frameId}.estSNR(target_idx)]) / length(target_idx);
    if LOG_ON
        disp(strcat("======>>>>> FrameId = ", num2str(frameId)));
        disp(strcat("======>>>>> Target Range bin = ", num2str(targetbin)));
        disp(strcat("======>>>>> Target Index Range bin = ", num2str(target_idx)));
        disp(strcat("======>>>>> Sub_intensity_SNR(dB) = ", num2str(Intensity_SNR(frameId))));
    end
    
    % Concatenate the phase data of each frame
    if (targetbin == 0)
        continue;
    else
        gPhase_RXs = cat(1, gPhase_RXs, squeeze(gRange_Profile(frameId, targetbin, :, :) ) );
    end
end
% Calculate the intensity signal-to-noise ratio
Intensity_SINR = sum(Intensity_SNR) / length(Intensity_SNR);
if LOG_ON
    disp(strcat("======>>>>> Intensity_SNR(dB) = ", num2str(Intensity_SNR)));
end

% IQ plan analysis the phase data
if LOG_ON
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
gCalib_Angle = frame_phase_calibration(gPhase_RXs, length(gframe_obj), NumRx, num_chirps, PLOT_ON);

% All RX antennas are averaged and the calibrated data signal is plotted.
gAngle = mean(gCalib_Angle(:, :), 2);
gAngle_diff = diff(gAngle);
if PLOT_ON
    figure(131)
    set(gcf,'units','normalized','outerposition',[0.1 0.4 0.6 0.4])
    subplot(121)
    plot(gAngle,'color', 'blue');
    title({"Phase Signal", strcat("Target Bin",num2str(targetbin))})
    subplot(122)
    plot(gAngle_diff,'color', 'blue');
    title("gBeams Phase Angle Diff")
    pause(0.01)
end

% Calculate second part Sensing Signal Intensity SNR
Phase_SINR = 0;




est_PSINR = Intensity_SINR + Phase_SINR;

end

