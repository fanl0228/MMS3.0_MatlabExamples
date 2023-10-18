function Hd = ellip_HPF
%ELLIP_HPF 返回离散时间滤波器对象。

% MATLAB Code
% Generated by MATLAB(R) 9.14 and DSP System Toolbox 9.16.
% Generated on: 14-Oct-2023 20:51:57

% Elliptic Highpass filter designed using FDESIGN.HIGHPASS.

% All frequency values are in Hz.
Fs = 2000;  % Sampling Frequency

Fstop = 0.1;     % Stopband Frequency
Fpass = 0.5;     % Passband Frequency
Astop = 80;      % Stopband Attenuation (dB)
Apass = 1;       % Passband Ripple (dB)
match = 'both';  % Band to match exactly

% Construct an FDESIGN object and call its ELLIP method.
h  = fdesign.highpass(Fstop, Fpass, Astop, Apass, Fs);
Hd = design(h, 'ellip', 'MatchExactly', match);

% [EOF]