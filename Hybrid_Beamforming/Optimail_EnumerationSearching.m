

tic
clearvars
close all
setenv('CASCADE_SIGNAL_PROCESSING_CHAIN_MIMO', 'D:\ti\mmwave_studio_03_00_00_14\mmWaveStudio\MMS3.0_MatlabExamples\Hybrid_Beamforming')

% Init Parameters of function Each_Steering_Calculate_pSINR 
dataFolder_Path = 'I:\20231014_EXP_LOS\LOS_Scence4_Sample02\'; 
isStepAngle = 0;
TxBF_Angle = 0;
RxBF_Angle = TxBF_Angle;
UsedFrames = 20;
NumRx = 16;
num_chirps = 252;
Signal_FS = 1000;
LOG_ON = 1;
PLOT_ON = 1;
png_floder = strcat([dataFolder_Path(1:end-1), '_png\']);
if ~exist(png_floder,'dir')
    mkdir(png_floder)
end
LogFileId = fopen(strcat(png_floder, 'LogFile.txt'), 'w');


rangeResolution = 0;   % calculate value
est_PSINR_all =[];
Intensity_estSINR_all =[];
Phase_estSINR_all = [];
range_angle_spec = [];

step = 1;
for TxBF_Angle = -60:60
    %
    
    fprintf(LogFileId, '==================================>>>>%s \t', dataFolder_Path);
    fprintf(LogFileId, '%s: %s, \n', "Beam steering angle (deg)", num2str(TxBF_Angle));
    
    %
    RxBF_Angle = TxBF_Angle;
    [Intensity_estSINR, Phase_estSINR, gframe_obj] ...
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
                                            LogFileId);
    
    
    Intensity_estSINR_all(step) = Intensity_estSINR;
    Phase_estSINR_all(step) = Phase_estSINR;
    est_PSINR_all(step) = Intensity_estSINR + Phase_estSINR;
    
    if ~isempty(gframe_obj)
        ra_spec = abs(gframe_obj{1}.RA_Spectrum);
        for j = 2:length(gframe_obj)
            ra_spec = ra_spec + abs(gframe_obj{j}.RA_Spectrum);
        end
        ra_spec = ra_spec/length(gframe_obj);
        ra_spec_log = 10*log10(ra_spec);
        range_angle_spec(:, step) = ra_spec_log;
    else
        range_angle_spec(:, step) = zeros(1, 128); % range bin size = 128
    end

    rangeResolution = gframe_obj{1}.rangeResolution;

    step = step +1;
    toc
end

%% -----------------------PSINR
fig225=figure(219);
search_angle = -floor(length(Intensity_estSINR_all)/2):1:floor(length(Intensity_estSINR_all)/2);
plot(search_angle, Intensity_estSINR_all, 'r')
hold on;
plot(search_angle, Phase_estSINR_all, 'b')
hold on;
plot(search_angle, est_PSINR_all, 'g')
hold off;
% calculate the beamangle
[val, idx] = max(est_PSINR_all);
BeamAngle = idx - 61;
%
xlim([-floor(length(Intensity_estSINR_all)/2), floor(length(Intensity_estSINR_all)/2)])
legend({"Intensity SINR", "Motion SINR", "PSINR"}, 'Location','best');
xlabel("Beamsearching Angle")
ylabel("SINR(dB)")
title({"Evaluate PSINR Metric", strcat("optimal BeamAngle: ", num2str(BeamAngle))})
pause(0.01)
% Save fig

psinr_file = strcat([png_floder, 'PSINR.png']);
saveas(fig225, psinr_file, 'png');
psinr_filefig = strcat([png_floder, 'PSINR.fig']);
saveas(fig225, psinr_filefig, 'fig');

fprintf(LogFileId, '================================================>>\n');
fprintf(LogFileId, '=Optimail TxBeamforming Angle (deg): ====>>>>>>%s deg: \n', num2str(BeamAngle));


fig220=figure(220);
PSINR2_phase = Phase_estSINR_all;
plot(search_angle, PSINR2_phase, 'k')
pause(0.05)

PSINR2_phase_file = strcat([png_floder, 'PSINR2_phase_file.png']);
saveas(fig220, PSINR2_phase_file, 'png');
psinr_filefig = strcat([png_floder, 'PSINR2_phase_file.fig']);
saveas(fig220, psinr_filefig, 'fig');

%% ----------------------极坐标显示
fig221=figure(221);
offset = 30;
theta = (offset+1:offset+length(Intensity_estSINR_all))/180*pi;
polarplot(theta, est_PSINR_all, 'b-','linewidth',1);
% pax = gca;
% pax.ThetaDir = 'clockwise';		   % 按顺时针方式递增
% pax.ThetaZeroLocation = 'top';     % 将0度放在顶部  
title("Beam pattern (dB)")
pause(0.01)
% Save fig
psinr_polar_floder = strcat([dataFolder_Path(1:end-1), '_png\']);
psinr_polar_file = strcat([psinr_polar_floder, 'psinr_polar.png']);
saveas(fig221, psinr_polar_file, 'png');
psinr_polar_filefig = strcat([psinr_polar_floder, 'psinr_polar.fig']);
saveas(fig221, psinr_polar_filefig, 'fig');


%% ---------------------- RxBF 可视化 RA 图
fig22=figure(222);
y_axis = (1:size(range_angle_spec,2))-61;
x_axis = ((1:size(range_angle_spec,1))-1)*rangeResolution;
fig222 = surf(y_axis, x_axis, range_angle_spec);
fig222.FaceColor = 'interp';
fig222.LineStyle="none";
colormap(gca,"jet");
title("RxBF Range Angle Heapmat")
xlabel("Agnel(deg)")
ylabel("Range(m)")
pause(0.01)

% Save fig
png_floder = strcat([dataFolder_Path(1:end-1), '_png\']);
rxbf_file = strcat([png_floder, 'RXBF.png']);
saveas(fig22, rxbf_file, 'png');
rxbf_filefig = strcat([png_floder, 'RXBF.fig']);
saveas(fig22, rxbf_filefig, 'fig');


fclose(LogFileId);
fclose all;
