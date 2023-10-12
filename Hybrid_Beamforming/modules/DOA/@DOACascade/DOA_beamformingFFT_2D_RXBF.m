
%DOA_beamformingFFT_2D_RXBF.m
%
%DOA_beamformingFFT_2D_RXBF function perform 2D angle estimation based on FFT beamforming, 
% both azimuth/elevation peak selection is done in 2D FFT domain


%input:
%   obj: object instance
%   sig: complex signal vector, with each value corresponding to each
%   antenna. The length of this vector equals to numTX x numRX enabled.
%   There can be overlapped antennas. this signal needs to be re-arranged
%   based on D value to form the virtual antenna array

%output:
%   angleObj_est: angle estimation results
%   angle_sepc_2D_fft: angle 2D fft spectrum

function [angleObj_est, angle_sepc_2D_fft, range_beam_angle]= DOA_beamformingFFT_2D_RXBF(obj, sig, current_obj, Txbeam_angle, DopplerFFTIn)

range_beam_angle = [];

%field of view to do beamforming
angles_DOA_az = obj.angles_DOA_az;
angles_DOA_ele = obj.angles_DOA_ele;

%distance unit in terms of wavelength
d = obj.antDis;
%2D matrix providing antenna coordinates

f0 = 77e9;
D_BF = obj.D(:,1);
angleFFTSize = obj.DOAFFTSize; 
DopplerFFTSize = obj.dopplerFFTSize;


%% ---test RxBF ----
if 1
    %decide non-zero doppler bins to be used for dynamic range-azimuth heatmap
    ratio = 0.5;
    DopplerPower = sum(mean((abs(DopplerFFTIn(:,:,:))),3),1);
    DopplerPower_noDC = DopplerPower([1: DopplerFFTSize/2-1 DopplerFFTSize/2+3:end]);
    [peakVal peakInd] = max(DopplerPower_noDC);
    threshold = peakVal*ratio;
    indSel = find(DopplerPower_noDC >threshold);
    for ii = 1:length(indSel)
        if indSel(ii) > DopplerFFTSize/2-1
            indSel(ii) = indSel(ii) + 3;
        end
    end
    % dynamic
    angle_range_dynamic = squeeze(sum(DopplerFFTIn(:, indSel, :), 2) );
    wx = sind(Txbeam_angle);
    a1_az = exp(1j*2*pi*d*(D_BF*wx));
    
    for irange = 1:size(angle_range_dynamic,1)
        RX_data = squeeze(angle_range_dynamic(irange, :));
        range_beam_angle(irange) = a1_az'*(RX_data'*RX_data)*a1_az;
    end
    % CFAR Detection
    [N_obj, Ind_obj, noise_obj, CFAR_SNR] = CFAR_SO(10*log10(abs(range_beam_angle).^2), 10, 5, 1.2, 0);
    if 0
        figure(333)
        plot(10*log10(abs(range_beam_angle).^2))
        hold on;
        scatter(Ind_obj, noise_obj)
        title("Range angle CFAR\_SO, N\_obj: ", num2str(N_obj));
        pause(0.01)
        hold off;
    end
end

%%
%FFT based implementation
%first form a 2D matrix based on the antenna coordinates
D = obj.D;
D = D + 1;
apertureLen_azim = max(D(:,1));
apertureLen_elev = max(D(:,2));
sig_2D = zeros(apertureLen_azim, apertureLen_elev);
for i_line = 1:apertureLen_elev
    ind = find(D(:,2) == i_line);
    D_sel = D(ind,1);
    sig_sel = sig(ind);
    [val, indU] = unique(D_sel);
    
    sig_2D(D_sel(indU), i_line) = sig_sel(indU);
    
    % RxBF
    wx = sind(Txbeam_angle);
    a1_az = exp(-1j*2*pi*f0*d*(D_BF*wx)); 
    sig_2D_RxBF(D_sel(indU), i_line) = sig_sel(indU) .* a1_az;  
      
end
% 
%%
% run FFT on azimuth and elevation
angle_sepc_1D_fft_RxBF = fftshift(fft(sig_2D_RxBF, angleFFTSize,1),1); 
angle_sepc_1D_fft_RxBF = flipud(angle_sepc_1D_fft_RxBF);
% angle_sepc_1D_fft = angle_sepc_1D_fft_RxBF;
% No RXBF
angle_sepc_1D_fft = fftshift(fft(sig_2D, angleFFTSize,1),1); 
angle_sepc_1D_fft = flipud(angle_sepc_1D_fft);

angle_sepc_2D_fft=fftshift(fft(angle_sepc_1D_fft, angleFFTSize, 2),2);  % 

if 0
    figure(565)
    plot(abs(angle_sepc_1D_fft),'b');
    hold on;
    plot(abs(angle_sepc_1D_fft_RxBF), 'r');
    hold off;
    pause(0.01)
end


wx_vec=[-pi:2*pi/angleFFTSize:pi];
wz_vec=[-pi:2*pi/angleFFTSize:pi];
wx_vec = wx_vec(1:end-1);
wz_vec = wz_vec(1:end-1);
%use one row with complete azimuth antenna of 1D FFT output for azimuth
%estimation
spec_azim = abs(angle_sepc_1D_fft_RxBF(:,1));
obj.sidelobeLevel_dB = obj.sidelobeLevel_dB_azim;

[peakVal_azim, peakLoc_azim] = DOA_BF_PeakDet_loc(obj, spec_azim);

if apertureLen_elev ==1
    %azimuth array only, no elevation antennas
    obj_cnt = 1;
    angleObj_est= [];
    for i_obj = 1:length(peakLoc_azim)
        ind = peakLoc_azim(i_obj);
        
        azim_est = asind(wx_vec(ind)/(2*pi*d)); % 估算目标所在的角度值
        if (azim_est >= angles_DOA_az(1) && azim_est <= angles_DOA_az(2))
            angleObj_est(1, obj_cnt) = azim_est;
            angleObj_est(2, obj_cnt) = 0;
            angleObj_est(3, obj_cnt) = ind;
            angleObj_est(4, obj_cnt) = 0;

            locs = setdiff(1:length(spec_azim), ind);
            angleObj_est(5, obj_cnt) = 10*log10((peakVal_azim(i_obj).^2)./mean(spec_azim(locs)));

            obj_cnt = obj_cnt+1;
            
        else
            continue;
        end
    end
    
else
    %azimuth and elevation angle estimation
   
    % figure(1);plot(spec_azim); hold on; grid on
    % plot(peakLoc_azim, spec_azim(peakLoc_azim),'ro');hold on
    
    %for each detected azimuth, estimate its elevation
    % figure(2)
    obj_cnt = 1;
    angleObj_est= [];
    obj.sidelobeLevel_dB = obj.sidelobeLevel_dB_elev;
    for i_obj = 1:length(peakLoc_azim)
        ind = peakLoc_azim(i_obj);
        spec_elev = abs(angle_sepc_2D_fft(ind,:));
        [peakVal_elev, peakLoc_elev] = DOA_BF_PeakDet_loc(obj, spec_elev);
        %calcualte the angle values
        for j_elev = 1:length(peakVal_elev)
            azim_est = asind(wx_vec(ind)/(2*pi*d));
            elev_est = asind(wz_vec(peakLoc_elev(j_elev))/(2*pi*d));
            
            if (azim_est >= angles_DOA_az(1) && azim_est <= angles_DOA_az(2) ...
                    &&elev_est >= angles_DOA_ele(1) && elev_est <= angles_DOA_ele(2))
                angleObj_est(1,obj_cnt) = azim_est;
                angleObj_est(2,obj_cnt) = elev_est;              
                angleObj_est(3,obj_cnt) = ind;
                angleObj_est(4,obj_cnt) = peakLoc_elev(j_elev);
                %plot(angleObj_est(4,obj_cnt),angleObj_est(3,obj_cnt) ,'x','MarkerSize',12, 'LineWidth',2);
                %hold on
                obj_cnt = obj_cnt+1;
                
            else
                continue;
            end
        end        
    end    
    %hold off

end










