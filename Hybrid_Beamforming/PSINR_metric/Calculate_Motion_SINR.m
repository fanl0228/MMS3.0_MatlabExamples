

% phaseAngle: [chirps*frames, nRx]
%
function [motion_estSNR] = Calculate_Motion_SINR(phaseAngle, Signal_FS, nFrames, nChirps, nRx)

% input assert
[cur_nSig, cur_nRx] = size(phaseAngle);
if (cur_nSig ~= (nChirps * nFrames)) || (cur_nRx ~= nRx)
    disp(["PhaseAngle data size is ERROR..."]);
    return;
end

% 
% All RX antennas are averaged and the calibrated data signal is plotted.
gAngle = mean(phaseAngle(:, :), 2);
gAngle_diff = diff(gAngle);

% 高通滤波去直流
gAngle_diff = filter(ellip_HPF, gAngle_diff);

% FFT
nFFTSize = 1024;
gAngle_diff_fft = abs(fftshift(fft(gAngle_diff, nFFTSize)));

% 频域直接去直流
gAngle_diff_fft(length(gAngle_diff_fft)/2+1) = 0;

gAngle_diff_fft_half = gAngle_diff_fft(nFFTSize/2+1:end);
gAngle_diff_fft_norm_half = gAngle_diff_fft_half / (nFFTSize/2);

if 0
    figure(154)
    plot(gAngle_diff_fft)
    pause(0.01)

    figure(158)
    plot(gAngle_diff_fft - mean(gAngle_diff_fft))
    pause(0.01)

    figure(159)
    plot(gAngle_diff_fft_norm_half)
    pause(0.01)
end

% find peaks
[val, idx] = max(gAngle_diff_fft_norm_half);
[pks, locs] = findpeaks(gAngle_diff_fft_norm_half, "MinPeakHeight", val*0.5);

freqs = (locs-1)*(Signal_FS/nFFTSize);

signal_power = sum(gAngle_diff_fft_norm_half(locs).^2) / length(locs);

% Power
noise_idx = setdiff( (1:(nFFTSize/2)),  [locs]);

noise_power = sum(gAngle_diff_fft_norm_half(noise_idx).^2) /length(noise_idx);

motion_estSNR = 10*log10(signal_power/noise_power);

%motion_estSNR = 10*log10(signal_power);

end

