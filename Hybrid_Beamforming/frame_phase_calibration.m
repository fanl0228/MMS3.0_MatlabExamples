function [gCalib_Angle] = frame_phase_calibration(gPhase_RXs, validFramsNum, numRxToEnable, nchirp_loops, PLOT_ON)
    
    if (~isempty(gPhase_RXs))

        for idxRX = 1:numRxToEnable
            % 1. get angle 
            Phase_Angle = unwrap(angle(gPhase_RXs(:, idxRX)));

            if 0    
                figure(7);
                plot(Phase_Angle);
                title("unwrap angle value")
                pause(0.01)
                
                signal = diff(Phase_Angle);
                signal = signal(1:length(signal)-1);

                figure(8);
                plot(signal);
                title("[Comparing] original phase diff signal,, No Smooth")
                pause(0.01)
            end
            
            % 2.fram 之间的 gap 相位不连续优化处理 
            if (validFramsNum> 2)
                % step1 去除 由 number chirp 引起的 相位突变
                angle_gap = 0.0;
                for idxChirp = 1:size(Phase_Angle, 1)/nchirp_loops-1 

                    angle_gap = angle_gap ...
                                + Phase_Angle(idxChirp * nchirp_loops + 1) ...
                                - Phase_Angle(idxChirp * nchirp_loops);
                    
                    Phase_Angle(idxChirp * nchirp_loops : ...
                                (idxChirp+1) * nchirp_loops - 1) ...
                                    = Phase_Angle(idxChirp * nchirp_loops : ...
                                                (idxChirp + 1) * nchirp_loops - 1) - angle_gap;
                end
                

                % step2
                for idxChirp = 1:size(Phase_Angle, 1)/nchirp_loops-1 
                    th1 = abs( Phase_Angle(idxChirp*nchirp_loops-1) ...
                                -  Phase_Angle(idxChirp*nchirp_loops));
                    th2 = abs( Phase_Angle(idxChirp*nchirp_loops+1) ...
                                -  Phase_Angle(idxChirp*nchirp_loops));
                    
                    Threshold = abs(Phase_Angle(idxChirp*nchirp_loops+1) ...
                                    - Phase_Angle(idxChirp*nchirp_loops-1));
            
                    if th1 + th2 > Threshold
                        Phase_Angle(idxChirp*nchirp_loops) ...
                                = 0.5*(Phase_Angle(idxChirp*nchirp_loops+1) ...
                                    + Phase_Angle(idxChirp*nchirp_loops-1));
                    end
                end
            
            end
            
            % 保存每一个RX 通道的相位数据，【考虑后续在RX通道做beamforming】
            gCalib_Angle(1:length(Phase_Angle), idxRX) = Phase_Angle;
            
            if false
                % step 相位差分
                Beam_Sgnal = diff(Phase_Angle);
                Beam_Sgnal = Beam_Sgnal(1:length(Beam_Sgnal) - 1);
                figure(111)
                subplot(121)
                plot(unwrap(angle(gPhase_RXs(:, idxRX))),'color', 'blue' );
                hold on;
                plot(Phase_Angle,'color', 'red');
                hold off;
                title({"Phase Signal"})
                subplot(122)
                plot(Beam_Sgnal);
                title("gBeams Phase Angle Diff")
                pause(0.01)
            end %if(DEBUG_PLOTS)
                
        end % idxRx    
    end %if (~isempty(gPhase_RXs))



    





end

