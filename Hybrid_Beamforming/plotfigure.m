close all;

%% PSINR
filename = 'I:\HyBF_Datasets\00_mmWaveRL_Datasets\train\NLOS12_DRLSample0.h5';
h5disp(filename)
Intensity_ = h5read(filename, '/rewards_Intensity');
Phase_ = h5read(filename, '/rewards_Phase');
PSINR_ = h5read(filename, '/rewards_PSINR');

% 2D scatter
fig226=figure(1);
sizeMarker = 50;
colorMarker = linspace(1,10,length(Intensity_)); 
scatter(Intensity_, Phase_, sizeMarker, colorMarker, 'filled');
hold on;
Intensity_max = max(Intensity_);
Intensity_min = min(Intensity_);
Phase_max = max(Phase_);
Phase_min = min(Phase_);
plot([Intensity_max*0.8, Intensity_max*0.8],[Phase_min*0.8, Phase_max*1.2], 'r--',LineWidth=1)
xlim([Intensity_min*0.9, Intensity_max*1.1])
ylim([Phase_min*0.9, Phase_max*1.1])
xlabel("SINR_{cs} (dB)")
ylabel("SINR_{bs} (dB)")
pause(0.05)
% ============================


figure(2)
Xangel = 1:121;
plot(Xangel, 10.^(Phase_/10), LineWidth=1.5)
xlim([0, 121])
xticklabels([-60, -40,-20, 0, 20, 40, 60])
xlabel('TxBF Angle (deg)')
ylabel('Entorpy_{bs}(dB)')


%
fig1=figure(3);
Xangel = 1:121;
plot(Xangel, Phase_, LineWidth=1.5)
hold on;
plot(Xangel, 22*(Intensity_/(max(Intensity_)-min(Intensity_))), LineWidth=1.5)
xlim([0, 121])
xticklabels([-60, -40,-20, 0, 20, 40, 60])
xlabel('TxBF Angle (deg)')
ylim([0, 50])
ylabel('SINR(dB)')
legend(["SINR_{bs}", "SINR_{cs}"])


%% ----------------------极坐标显示
fig221=figure(4);
offset = 30;
theta = (offset+1:offset+length(Phase_))/180*pi;
polarplot(theta, Phase_,'linewidth',1);
thetalim([0,180]);     %修改此处，决定扇形区域
% pax = gca;
% pax.ThetaDir = 'clockwise';		   % 按顺时针方式递增
% pax.ThetaZeroLocation = 'top';     % 将0度放在顶部  
% title("Beam pattern (dB)")
pause(0.01)


%% RDA metric
State_RDA_Struct = load('H:\NLOS12_sample03\BeamAngle055_Sample03_State_RDA.mat');
State_RDA = State_RDA_Struct.State_RDA;
size(State_RDA)

% angle FFT
State_RDA_4dfft = fftshift(fft(State_RDA, 128, 4), 4);   % [frame, range, doppler, angle]
size(State_RDA_4dfft)

for VIS_FRAME = 1
    
    % figure Range-doppler
    RD_data = squeeze(mean(State_RDA_4dfft(VIS_FRAME, :, :, :), 4));
    figure(5)
    [y_axis, x_axis]=meshgrid(1: size(RD_data, 2), 1: size(RD_data, 1));
    fig2=surf(y_axis, x_axis, 10*log10(abs( RD_data ) ) );
    fig2.FaceColor = 'interp';
    fig2.LineStyle="none";
    colormap(gca,"jet");
    c = colorbar;
    c.Label.String = 'Relative Power(dB)';
    title("Range-Doppler")
    pause(0.05)
    
    
    % figure doppler-angle
    RDA_data = squeeze(State_RDA_4dfft(VIS_FRAME, :, :, :)); % r,d,a
    figure(6)
    [angle_axis, doppler_axis]=meshgrid(1: size(RDA_data, 3), 1: size(RDA_data, 2));
    fig1=surf(angle_axis, doppler_axis, 10*log10(abs( squeeze( RDA_data(1, :, :) ) ) ) );
    fig1.FaceColor = 'interp';
    fig1.LineStyle="none";
    colormap(gca,"jet");
    c = colorbar;
    c.Label.String = 'Relative Power(dB)';
    title("Range-Angle")
    pause(0.05)
    
    
    % range - angle - doppler
    % 创建一个随机的三维数组，你可以用你自己的数组代替
    figure(7)
    data = squeeze( 10*log10( abs( State_RDA_4dfft(VIS_FRAME, :, 33:96, 33:96) ) ) );
    RAD_data = permute(data,[1, 3, 2]);
    [X,Y,Z] = meshgrid(-31:1:32);
    
    x_slice = -1;   % 
    y_slice = -17;  %
    z_slice = [-8];
    
    Slice_fig=slice(X,Y,Z, RAD_data, x_slice, y_slice, z_slice,"nearest");
    
    Slice_fig(1).FaceColor = 'interp';
    Slice_fig(1).LineStyle="none";
    Slice_fig(2).FaceColor = 'interp';
    Slice_fig(2).LineStyle="none";
    Slice_fig(3).FaceColor = 'interp';
    Slice_fig(3).LineStyle="none";
    
    colormap(gca,"jet");
    c = colorbar;
    c.Label.String = 'Relative Power(dB)';
    
    xlabel("Angle")
    ylabel("Range")
    zlabel("Doppler")
    yticklabels([0,20,40,60,80])
    
    title("Slice of RDA")
    pause(0.05)
end





