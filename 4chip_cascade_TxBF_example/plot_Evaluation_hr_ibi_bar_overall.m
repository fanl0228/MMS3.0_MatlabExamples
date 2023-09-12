
%%
clc;
clf;
close all

%% HR
HR_mmWave = [63.75	64.25	63	64	63.25	63.25	64.25	63.25	64	63.25	62.75	63.75	64.25	65.25	65.5	64.75	62	60.75	61	58.75	60;
             90.75	91.25	94.5	86.25	91.25	92.75	92.5	94.5	92.25	92	89.75	89.5	91.75	91.5	89.25	94.25	94	96.25	96.25	98	97.5
             82.75	81	80.25	83.5	81.5	82.5	83.25	83	84.5	85.25	85.5	87	86.75	89.5	89.25	89.75	89.25	89.75	90.5	89	89.5;
             87	86	86	87	90	85.5	89.5	92.5	87	89	90.5	86	89	90	86.5	89	88	85.5	88.5	86.5	86.5]

HR_ECG = [64	64.25	63	64	63.5	63.5	64.5	63	64.5	63.5	63.25	64.25	64.25	65	65.25	65.25	61.75	60.25	61.5	59.25	60.5;
        91.25	90.75	93.75	85.75	90.5	92	92.75	93.75	92.75	92.25	89	89	92	91.75	89.5	94.5	93.5	96.75	95.75	98.5	98.5
        82.75	81	80	83	81.25	82.4	83	83	84.75	85.25	85.75	87.25	87	89.75	89.75	89.5	89.25	90	90.25	89.5	90;
        86.5	86	85.5	87.5	89.5	85	89.5	92	87	89.5	90	86.5	89.5	90	86	89.5	88.5	86	89	86.5	86 ]



angle = [-56.25	-50.625	-45	-39.375	-33.75	-28.125	-22.5	-16.875	-11.25	-5.625	0	5.625	11.25	16.875	22.5	28.125	33.75	39.375	45	50.625	56.25]

HR_Error_baseline = [1.5 1.7 0.9 1.1 1.2 0.6 0.8 0.4 0.7 0.6 0.5 0.7 0.6 0.9 0.7 0.6 1.2 1.3 1.1 1.5 1.3]
err0 = [0.05 0.07 0.06 0.09 0.07 0.06  0.09 0.01 0.02 0.06 0.08 0.04 0.07 0.06 0.08 0.06 0.08 0.1 0.12 0.1 0.09]


HR_Error = 100*abs(HR_mmWave - HR_ECG)./HR_ECG;
err1 = var(HR_Error);

HR_Error_baseline_mean = mean(HR_Error_baseline)
HR_Error_mean = mean(mean(HR_Error))

%% IBI
IBI_mmWave = [938.223	931.141	951.438	935.559	944.256	944.788	930.484	947.616	936.462	946.714	954.345	939.939	928.964	915.469	913.788	922.466	967.465	985.597	980.533	1019.495 997.318;
               659.056	655.505	633.358	595.226	656.884	646.395	648.125	633.912	649.162	651.681	667.436	669.78	652.835	654.807	670.973	636.439	636.827	622.159	682.657	666.198	654.039;
               723.442	739.265	744.441	717.169	734.711	725.71	720.345	722.073	709.588	701.934	700.95	687.078	689.941	669.291	672.138	666.652	671.096	667.853	662.258	673.025	670.068;
               672.665	681.541	682.922	686.656	665.469	699.123	668.226	644.916	685.56	669.979	660.808	695.712	673.525	665.399	690.21	672.454	676.687	696.564	663.047	681.087	681.033]

IBI_ECG = [931.969	932.358	945.81	933.781	943.685	942.394	926.233	951.619	926.729	941.953	945.328	927.969	929.276	919.292	917.67	919.326	969.911	996.249	970.537	1017.668	988.992;
           616	    624.125	605.873	585.739	634.909	589	625.717	613.662	657.077	596.229	592.297	595.98	605.691	600.742	608.102	603.708	617.412	608.041	600.06	607.533	606.741;
           743.492	758.123	768.425	726	739.446	728.848	723.518	721.735	707.894	700.622	696.047	683.347	691.287	670.128	669.081	670.503	673.524	666.889	667.634	673.944	670.533;
           686.705	690.279	692.959	683.611	666.324	699.906	663.642	651.391	683.402	665.922	663.422	690.127	665.922	664.444	690.884	665.922	673.763	696	667.236	686.705	690.605]
IBI_Error = 100*abs(IBI_ECG -IBI_mmWave)./IBI_ECG
% 计算 Ours 所有角度下的平均误差
IBI_Error_mean = mean(mean(IBI_Error))

% baseline
IBI_Error_baseline = 100*[53 67 59 41 52 36 28 34 27 36 35 27 36 29 27 26 37 49 58 55 64] ./ mean(IBI_ECG);
IBI_Error_baseline_mean = mean(IBI_Error_baseline)

data = [HR_Error_baseline_mean, HR_Error_mean; IBI_Error_baseline_mean, IBI_Error_mean]

figure('color','w')
color1=[46, 114, 188]/255;
color2=[206, 85, 30]/255;

h = bar(data)

for i = 1:2
    text(i-0.15, data(i, 1), num2str(data(i,1),'%.2f'),...
        'HorizontalAlignment','center',...
        'VerticalAlignment','bottom',...
        'FontSize',16, 'FontName','Times New Roman')
    text(i+0.15, data(i, 2), num2str(data(i,2), '%.2f'),...
        'HorizontalAlignment','center',...
        'VerticalAlignment','bottom',...
        'FontSize',16, 'FontName','Times New Roman')

end


set(h(1), 'FaceColor', color1)
set(h(2), 'FaceColor', color2)
set(gca, 'linewidth', 0.75);
set(gca, 'Fontname', 'Times New Roman', 'fontsize', 16);
set(get(gca, 'XLabel'), 'FontSize', 16);
set(get(gca, 'YLabel'), 'FontSize', 16);
set(get(gca, 'TITLE'), 'FontSize', 16);
set(gca, 'XTickLabel', {'HR', 'IBI'},'FontSize', 16, 'FontName', 'Times New Roman')
xlabel('', 'Fontname','Times New Roman', 'Fontsize',16)
ylabel('Estimation Error (%)', 'Fontname','Times New Roman', 'Fontsize',16)
ylim([0,8])
set(gca, 'ytick', 0:2:8)


legend({'Baseline', 'Ours'}, 'location', 'NorthEast', 'FontSize', 14)




