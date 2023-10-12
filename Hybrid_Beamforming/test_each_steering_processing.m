

tic
clearvars
close all
setenv('CASCADE_SIGNAL_PROCESSING_CHAIN_MIMO', 'D:\ti\mmwave_studio_03_00_00_14\mmWaveStudio\MMS3.0_MatlabExamples\Hybrid_Beamforming')

UsedFrames = 10;
NumRx = 16;
num_chirps = 252;
LOG_ON = 1;
PLOT_ON = 0;


dataFolder_Path = 'I:\20230928_Exp_NLOS\Indoor_Scene4\test1\'; 

i = 1;
est_PSINR_all =[];
Intensity_estSINR_all =[];
Phase_estSINR_all = [];

range_angle_spec = [];

for Txbeam_angle = -60:60

    [Intensity_estSINR, Phase_estSINR, gframe_obj] = Each_Steering_Calculate_pSINR(dataFolder_Path, Txbeam_angle, UsedFrames, NumRx, num_chirps, LOG_ON, PLOT_ON);
    Intensity_estSINR_all(i) = Intensity_estSINR;
    Phase_estSINR_all(i) = Phase_estSINR;
    est_PSINR_all(i) = Intensity_estSINR + Phase_estSINR;
    
    ra_spec = abs(gframe_obj{1}.RA_Spectrum);
    for j = 2:length(gframe_obj)
        ra_spec = ra_spec + abs(gframe_obj{j}.RA_Spectrum);
    end
    ra_spec = ra_spec/length(gframe_obj);
    ra_spec_log = 10*log10(ra_spec);
    range_angle_spec(:, i) = ra_spec_log;

    i = i +1;
    toc
end

% -----------------------PSINR
figure(225)
search_angle = -floor(length(Intensity_estSINR_all)/2):1:floor(length(Intensity_estSINR_all)/2);
plot(search_angle, flip(Intensity_estSINR_all), 'r')
hold on;
plot(search_angle, flip(Phase_estSINR_all), 'b')
hold on;
plot(search_angle, flip(est_PSINR_all), 'g')
% hold on;
% PSINR2_power = 10*log10(10.^(Intensity_estSINR_all./10) + 10.^(Phase_estSINR_all./10));
% plot(search_angle, flip(PSINR2_power), 'k')
hold off;
xlim([-floor(length(Intensity_estSINR_all)/2), floor(length(Intensity_estSINR_all)/2)])
legend(["Intensity SINR", "Motion SINR", "PSINR1", "PSINR2"]);
xlabel("Beamsearching Angle")
ylabel("SINR(dB)")
title("Evaluate PSINR Metric")
pause(0.01)

% ----------------------极坐标显示
figure(221)
offset = 30;
theta = (offset+1:offset+length(Intensity_estSINR_all))/180*pi;
%polarplot(theta, flip(est_PSINR_all/max(est_PSINR_all)), 'b-','linewidth',1);
polarplot(theta, flip(est_PSINR_all), 'b-','linewidth',1);
% pax = gca;
% pax.ThetaDir = 'clockwise';		   % 按顺时针方式递增
% pax.ThetaZeroLocation = 'top';     % 将0度放在顶部  
title("Beam pattern (dB)")
pause(0.01)

% ---------------------- RxBF 可视化 RA 图
figure(222)
y_axis = 1:size(range_angle_spec,2);
x_axis = 1:size(range_angle_spec,1);
fig222 = surf(y_axis, x_axis, range_angle_spec);
fig222.FaceColor = 'interp';
fig222.LineStyle="none";
colormap(gca,"jet");
pause(0.01)



Intensity_estSINR_all1 = 10.^(Intensity_estSINR_all./10);
Intensity_estSINR_all11 = Intensity_estSINR_all1./max(Intensity_estSINR_all1);

SINR = Intensity_estSINR_all11 .* Phase_estSINR;



