close all;
%% Performance of PSINR Useability
PSINR_all_samples_Xdata = [];
PSINR_all_samples_Ydata = [];

Phase_all_samples_Xdata = [];
Phase_all_samples_Ydata = [];

Intensity_all_samples_Xdata = [];
Intensity_all_samples_Ydata = [];
for i = 1:5
    path = strcat(['I:\HyBF_Datasets\NLOS21_DRL\LOS21_Sample', num2str(i),'_PSINR.fig']);
    open(path);
    
    sample_all_lines = findall(gca, 'type', 'line');                    % 如果图中有多条曲线，lh为一个数组
    sample_all_lines_Xdata = get(sample_all_lines, 'xdata');            % 取出x轴数据，xc是一个元胞数组
    sample_all_lines_Ydata = get(sample_all_lines, 'ydata');            % 取出y轴数据，yc是一个元胞数组
    
    % Use line data
    PSINR_all_samples_Xdata(i,:) = sample_all_lines_Xdata{1};
    PSINR_all_samples_Ydata(i,:) = sample_all_lines_Ydata{1};%*35/max(sample_all_lines_Ydata{1});

    Phase_all_samples_Xdata(i,:) = sample_all_lines_Xdata{2};
    Phase_all_samples_Ydata(i,:) = sample_all_lines_Ydata{2};

    Intensity_all_samples_Xdata(i,:) = sample_all_lines_Xdata{3};
    Intensity_all_samples_Ydata(i,:) = sample_all_lines_Ydata{3};
end

figure()
search_angle = -floor(size(PSINR_all_samples_Ydata, 2)/2):1:floor(size(PSINR_all_samples_Ydata, 2)/2);
plot(PSINR_all_samples_Xdata', PSINR_all_samples_Ydata', LineWidth=1);
hold on;
plot(search_angle, mean(PSINR_all_samples_Ydata), LineWidth=1);
hold on;
envelope = max(PSINR_all_samples_Ydata);
%envelope = medfilt1(envelope, 5);
span = 0.1;                             % 适当调整span的值
envelope = smooth(envelope, span, 'loess');
plot(search_angle, envelope, 'r-', LineWidth=2);
xlim([-floor(size(PSINR_all_samples_Ydata, 2)/2), floor(size(PSINR_all_samples_Ydata, 2)/2)])
%ylim([0,42])
xlabel("Beamsearching Angle")
ylabel("PSINR (dB)")
xticks(-60:20:60)
yticklabels([-10, 0, 10, 20, 30, 40, 50])
grid on;
%title({"Evaluate PSINR Metric", strcat("optimal BeamAngle: ", num2str(BeamAngle))})
pause(0.01)


%% 

color=[
92 158 173
210 204 161 
206 190 190
237 177 131
239 111 108
];

figure();
sizeMarker = 30;
for i = 1:5
    scatter(Intensity_all_samples_Ydata(i,:), Phase_all_samples_Ydata(i, :), sizeMarker, 'filled');
    hold on;
end



