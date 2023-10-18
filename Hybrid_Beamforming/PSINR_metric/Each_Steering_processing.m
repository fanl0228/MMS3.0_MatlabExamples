% ***** output
% gframe_obj:
% gRange_Profile: [frame, range, chirps, rx]
% gRangeDoppler_Profile:


function [gframe_obj] = Each_Steering_processing(dataFolderName, TxBF_Angle, RxBF_Angle, UsedFrames, LOG_ON, PLOT_ON)
    %tic

    dataPlatform = 'TDA2';
    pro_path = getenv('CASCADE_SIGNAL_PROCESSING_CHAIN_MIMO');
    input_path = strcat(pro_path,'\utils\cascade_json_parser\');
    pathGenParaFile = [input_path,'generate_HybridBF_param.m'];
    clear(pathGenParaFile); % important to clear the same.m file, since Matlab does not clear cache automatically
    module_param_file = [input_path, 'module_param.m']; 
    
    gframe_obj = {};    


    % 判断文件夹是否存在
    [fileIdx_unique] = getUniqueFileIdx(dataFolderName);
        
    [fileNameStruct] = getBinFileNames_withIdx(dataFolderName, fileIdx_unique{1});

    [numValidFrames, ~] = getValidNumFrames(fullfile(dataFolderName, fileNameStruct.masterIdxFile));
    if numValidFrames < UsedFrames
        UsedFrames = numValidFrames;
        disp(strcat("UsedFrames is Changed >>>>> ", num2str(UsedFrames)));  
    end

    %generate parameter file for the dataset
    parameter_file_gen_AngleSweep_json(dataFolderName, module_param_file, pathGenParaFile, dataPlatform);

    %generate calibration matrix object for the dataset
    genCalibrationMatrixObj      = genCalibrationMatrixCascade('pfile', pathGenParaFile, ...
                                                                'calibrateFileName', dataFolderName);
    genCalibrationMatrixObj.binDataFile = fileNameStruct;

    % simTopObj is used for top level parameter parsing and data loading and saving
    simTopObj                   = simTopCascade('pfile', pathGenParaFile);
    rangeFFTObj                 = rangeProcCascade('pfile', pathGenParaFile);
    DopplerFFTObj               = DopplerProcClutterRemove('pfile', pathGenParaFile);
    CFAR_DetectionObj           = CFAR_CASO('pfile', pathGenParaFile);
    DOAObj                      = DOACascade('pfile', pathGenParaFile);

    
    valid_obj_frameId = 0;
    % get frame 
    for frameId = 1: UsedFrames
        genCalibrationMatrixObj.frameIdx = frameId;
        % force cal object to single TX phase/frame (numChirpPerLoop <-> numTXPhases) 
        genCalibrationMatrixObj.TxToEnable = 1;
        % read in .bin files， rawData = [samples, chirps, Rx]
        rawData = cascade_Read_rawData(genCalibrationMatrixObj);

        %% -----------------Range FFT------------------------------
        %perform 2D FFT
        rangeFFTOut = [];
        DopplerFFTOut = [];
        
        % range FFT
        rangeFFTOut(:,:,:)     = datapath(rangeFFTObj, rawData(:,:,:));
        % visualization Range FFT Output
        if false
            figure(100)
            plot(10*log10(abs(rangeFFTOut(:,:,1))));
            title("Plot Range FFT Output");
            pause(0.01);
        end

        
        %% -----------------Doppler FFT------------------------------
        DopplerFFTOut(:,:,:)   = datapath(DopplerFFTObj, rangeFFTOut(:,:,:));

        
        % Visualization Doppler FFT Output
        if false
            figure(101)
            fig2 = pcolor(10*log10( abs( DopplerFFTOut(:,:,1) ) ) );
            fig2.FaceColor = 'interp';
            fig2.LineStyle="none";
            title("Imagesc Range-Doppler FFT Output");
            pause(0.01)
        end

        % RD Signature log
        RD_Sig_integrate = 10*log10(sum((abs(DopplerFFTOut)).^2, 3) + 1);

        %% -----------------2D CFAR Detection------------------------------
        CFAR_detection_results = datapath(CFAR_DetectionObj, DopplerFFTOut);

        detect_all_points = [];
        for iobj = 1:length(CFAR_detection_results)
            detect_all_points(iobj, 1)      =CFAR_detection_results(iobj).rangeInd;
            detect_all_points(iobj, 2)      =CFAR_detection_results(iobj).dopplerInd_org;
            detect_all_points(iobj, 4)      =CFAR_detection_results(iobj).estSNR;
        end

        %% -----------------Range Angle, Point Cloud Generation------------
        angles_all_points = [];
        xyz = [];
        if ~isempty(CFAR_detection_results)
            % DOA, the results include detection results + angle estimation results.
            % access data with angleEst{frame}(objectIdx).fieldName
           
            % sort 
            [~, estSNR_Rank_Idxs] = sort([CFAR_detection_results.estSNR], "descend");
            if (length(estSNR_Rank_Idxs) > 10)
                SNR_TOPN = 10;
            else
                SNR_TOPN = length(estSNR_Rank_Idxs);
            end

            % dellet data of CFAR_detection_results.doppler=0
            NoDopplerZero_Idxs = find([CFAR_detection_results.doppler]); % find doppler non-zero

            estSNR_Rank_Idxs = intersect(estSNR_Rank_Idxs(1:SNR_TOPN), NoDopplerZero_Idxs, "stable");
                    
            if isempty(estSNR_Rank_Idxs)
                if LOG_ON
                    disp(["estSNR_Rank_Idxs is empty, Skip FrameId: ", num2str(frameId)])
                end
                continue;      
            end

            new_iobj = 0;
            for iobj = estSNR_Rank_Idxs
                new_iobj = new_iobj  + 1;
                CFAR_detection_results_rank(1, new_iobj).rangeInd                    = CFAR_detection_results(1, iobj).rangeInd;
                CFAR_detection_results_rank(1, new_iobj).range                       = CFAR_detection_results(1, iobj).range;
                CFAR_detection_results_rank(1, new_iobj).dopplerInd_org              = CFAR_detection_results(1, iobj).dopplerInd_org;
                CFAR_detection_results_rank(1, new_iobj).dopplerInd                  = CFAR_detection_results(1, iobj).dopplerInd;
                CFAR_detection_results_rank(1, new_iobj).doppler                     = CFAR_detection_results(1, iobj).doppler;
                CFAR_detection_results_rank(1, new_iobj).doppler_corr                = CFAR_detection_results(1, iobj).doppler_corr;
                CFAR_detection_results_rank(1, new_iobj).noise_var                   = CFAR_detection_results(1, iobj).noise_var;
                CFAR_detection_results_rank(1, new_iobj).bin_val                     = CFAR_detection_results(1, iobj).bin_val;
                CFAR_detection_results_rank(1, new_iobj).bin_val_org                 = CFAR_detection_results(1, iobj).bin_val_org;
                CFAR_detection_results_rank(1, new_iobj).estSNR                      = CFAR_detection_results(1, iobj).estSNR;
                CFAR_detection_results_rank(1, new_iobj).doppler_corr_overlap        = CFAR_detection_results(1, iobj).doppler_corr_overlap;
                CFAR_detection_results_rank(1, new_iobj).doppler_corr_FFT            = CFAR_detection_results(1, iobj).doppler_corr_FFT;
            end

            %% -----------------DOA Processing -----------------------------
            angleEst = datapath(DOAObj, CFAR_detection_results_rank, TxBF_Angle, RxBF_Angle, DopplerFFTOut);
            
            maxRangeShow = CFAR_DetectionObj.rangeBinSize * rangeFFTObj.rangeFFTSize;

            if ~isempty(angleEst)
                %% Filting using the doppler, estSNR, and Range angle fft
                % 去掉doppler_corr为0的目标,并且SNR取 top10
                idxs0 = find([angleEst.rangeInd] < rangeFFTObj.rangeFFTSize/2); % find range bin < rangeFFTSize/2

                idxs1 = find([angleEst.doppler_corr]); % find doppler non-zero
                
                idxs1 = intersect(idxs0, idxs1, "stable");
                
                if 0
                    Obj_rangeInds = [angleEst(1, idxs1).rangeInd];
                    
                    Angle_range = mag_data_dynamic(:, Obj_rangeInds, valid_obj_frameId); % [angle, range, frame]
                    %Angle_range(:,1:size(RA_data,2)) = RA_data(:,:);
                    if 0
                        figure(22)
                        fig2 = pcolor(Angle_range);
                        fig2.FaceColor = 'interp';
                        fig2.LineStyle="none";
                        colormap(gca,"jet");
                        c = colorbar;
                        c.Label.String = 'Relative Power(dB)';                    
                    end
                    
                    [vals, locs] = max(max(Angle_range));
                    Angle_range(Angle_range<0.8*vals) = 0;
    
                    if 0
                        figure(23)
                        fig2 = pcolor(Angle_range);
                        fig2.FaceColor = 'interp';
                        fig2.LineStyle="none";
                        colormap(gca,"jet");
                        c = colorbar;
                        c.Label.String = 'Relative Power(dB)';                    
                    end
    
                    [rows, cols]=find(Angle_range(:,:)~=0);
                end

                % [~, idxs2] = sort([angleEst.estSNR], "descend");
                % if (length(idxs2) > 20)
                %     SNR_TOPN = 20;
                % else
                %     SNR_TOPN = length(idxs2);
                % end
                % idxs3_array = intersect(idxs1, idxs2(1:SNR_TOPN));
                       
                %idxs3_array = unique(cols);
                idxs3_array = idxs1;

                if isempty(idxs3_array)
                    if LOG_ON
                        disp(["idxs3 is empty, and no object can be detected."])
                    end
                    continue;
                end
                
                % Use all cfar results
                % idxs3 = intersect(idxs1, idxs2);

                % 
                % 计算target信噪比，从中取SNR 前top N个点作为基础点。 angleEst
                new_iobj = 0;
                for iobj = idxs3_array
                    new_iobj = new_iobj  + 1;
                    out_frame(1, new_iobj).rangeInd                    = angleEst(1, iobj).rangeInd;
                    out_frame(1, new_iobj).dopplerInd                  = angleEst(1, iobj).dopplerInd;
                    out_frame(1, new_iobj).range                       = angleEst(1, iobj).range;
                    out_frame(1, new_iobj).doppler                     = angleEst(1, iobj).doppler;
                    out_frame(1, new_iobj).doppler_corr                = angleEst(1, iobj).doppler_corr;
                    out_frame(1, new_iobj).dopplerInd_org              = angleEst(1, iobj).dopplerInd_org;
                    out_frame(1, new_iobj).noise_var                   = angleEst(1, iobj).noise_var;
                    out_frame(1, new_iobj).bin_val                     = angleEst(1, iobj).bin_val;
                    out_frame(1, new_iobj).estSNR                      = angleEst(1, iobj).estSNR;
                    out_frame(1, new_iobj).doppler_corr_overlap        = angleEst(1, iobj).doppler_corr_overlap;
                    out_frame(1, new_iobj).doppler_corr_FFT            = angleEst(1, iobj).doppler_corr_FFT;
                    out_frame(1, new_iobj).angles                      = angleEst(1, iobj).angles;
                    out_frame(1, new_iobj).spectrum                    = angleEst(1, iobj).spectrum;
                    out_frame(1, new_iobj).range_beam_spectrum         = angleEst(1, iobj).range_beam_spectrum;
                end

                valid_obj_frameId = valid_obj_frameId + 1;
                detection_results_all_frames{valid_obj_frameId} = out_frame;
                
                % calculate target azimuth
                % [~, idx_angle] = max( 10*log10( abs( [angleEst.range_beam_spectrum] ) ) );
                % if(ismember(idx_angle, idxs3_array))
                %     obj_rangebin = idx_angle;
                % else
                %     obj_rangebin= [detection_results_all_frames{valid_obj_frameId}.rangeInd];
                % end
                obj_rangebin= [detection_results_all_frames{valid_obj_frameId}.rangeInd];

                % [Beam, frame, objs_range_bin]
                gframe_obj_rangeBins(valid_obj_frameId, 1:length(obj_rangebin)) = obj_rangebin; 
                gframe_obj_estSNR(valid_obj_frameId, 1:length(obj_rangebin))  = detection_results_all_frames{valid_obj_frameId}.estSNR;
                gframe_obj_estAngle(valid_obj_frameId, :, 1:length(obj_rangebin)) = [detection_results_all_frames{valid_obj_frameId}.angles];
                gframe_obj_RA_Spectrum(valid_obj_frameId, :)  = detection_results_all_frames{valid_obj_frameId}.range_beam_spectrum;
                
                            
                % Visualization
                if (PLOT_ON && (~mod(frameId,5)))
                    fig103=figure(103);
                    set(gcf,'units','normalized','outerposition',[0.05 0.3 0.6 0.6])                
                    
                    subplot(2,3,1)               
                    plot((1:size(RD_Sig_integrate,1))*CFAR_DetectionObj.rangeBinSize, RD_Sig_integrate(:,size(RD_Sig_integrate,2)/2+1),'g','LineWidth',4);
                    hold on; grid on;
                    for ii=1:size(RD_Sig_integrate,2)
                        plot((1:size(RD_Sig_integrate,1))*CFAR_DetectionObj.rangeBinSize, RD_Sig_integrate(:, ii));
                        hold on; grid on;
                        if ~isempty(CFAR_detection_results)
                            ind = find(detect_all_points(:,2) == ii);
                            if (~isempty(ind))


                                rangeInd = detect_all_points(ind,1);
                                plot(rangeInd*CFAR_DetectionObj.rangeBinSize, RD_Sig_integrate(rangeInd, ii), ...
                                        'o','LineWidth',2,...
                                        'MarkerEdgeColor','k',...
                                        'MarkerFaceColor',[.49 1 .63],...
                                        'MarkerSize',10);
                            end
                        end
                    end
                    
                    xlabel('Range(m)');
                    ylabel('Receive Power (dB)')
                    title({'Range Profile(zero Doppler - thick green line):',' valid obj frameId ' num2str(valid_obj_frameId)});
                    hold off;

                    subplot(2,3,2);
                    %subplot_tight(2,2,2,0.1)
                    RD_Sig_integrate_noDC = RD_Sig_integrate(:, [1:size(RD_Sig_integrate,2)/2-1 size(RD_Sig_integrate,2)/2+3:end]);
                    [Y, X]=meshgrid(1:1:size(RD_Sig_integrate_noDC,1), -size(RD_Sig_integrate_noDC,2)/2+1:1:size(RD_Sig_integrate_noDC,2)/2);
                    Y = Y * CFAR_DetectionObj.rangeBinSize;
                    X = X * CFAR_DetectionObj.velocityBinSize;
                    fig2 = pcolor(X,Y, RD_Sig_integrate_noDC');
                    fig2.FaceColor = 'interp';
                    fig2.LineStyle="none";
                    colormap(gca,"jet");
                    c = colorbar;
                    c.Label.String = 'Relative Power(dB)';
                
                    ylim([1 maxRangeShow/2])

                    [RD_SNR_max_val, RD_SNR_max_idx]= max(detect_all_points(:, 4));
                    title({' Range/Velocity Plot No-DC', ...
                            strcat('obj1 SNR=', num2str(CFAR_detection_results(1).estSNR)) ...
                            strcat(['SNR_{max}(dB)' num2str(RD_SNR_max_val) ...
                                    ' Range(m)=' num2str(detect_all_points(RD_SNR_max_idx, 1)*CFAR_DetectionObj.rangeBinSize)] )  ...
                            });
                    pause(0.01)
 
                    % plot 3D point cloud
                    if length(out_frame) > 0
                        for iobj = 1:length(out_frame)
                            angles_all_points(iobj,1:2)=out_frame(iobj).angles(1:2);
                            angles_all_points(iobj,3)=out_frame(iobj).estSNR;
                            angles_all_points(iobj,4)=out_frame(iobj).rangeInd;
                            angles_all_points(iobj,5)=out_frame(iobj).doppler_corr;
                            angles_all_points(iobj,6)=out_frame(iobj).range;
                            %switch left and right, the azimuth angle is not flipped
                            xyz(iobj,1) = angles_all_points(iobj,6) * sind( angles_all_points(iobj,1) * -1 ) * cosd( angles_all_points(iobj,2) );
                            xyz(iobj,2) = angles_all_points(iobj,6) * cosd( angles_all_points(iobj,1) * -1 ) * cosd( angles_all_points(iobj,2) );
                            %switch upside and down, the elevation angle is flipped
                            xyz(iobj,3) = angles_all_points(iobj,6) * sind(angles_all_points(iobj,2) * -1);
                            xyz(iobj,4) = out_frame(iobj).doppler_corr;
                            xyz(iobj,9) = out_frame(iobj).dopplerInd_org;
                            xyz(iobj,5) = out_frame(iobj).range;
                            xyz(iobj,6) = out_frame(iobj).estSNR;
                            xyz(iobj,7) = out_frame(iobj).doppler_corr_overlap;
                            xyz(iobj,8) = out_frame(iobj).doppler_corr_FFT; 
                        end

                        angles_all_all{valid_obj_frameId} = angles_all_points;
                        xyz_all{valid_obj_frameId}  = xyz;
                        
                        %tic
                        
                        if PLOT_ON
                            moveID = find(abs(xyz(:,4))>=0);
                            subplot(2,3,6);                        
                            
                            if valid_obj_frameId==1
                                scatter3(xyz(moveID,1),xyz(moveID,2),xyz(moveID,3),45,(xyz(moveID,4)), 'filled');
                                
                            else
                                yz = [xyz_all{valid_obj_frameId}; xyz_all{valid_obj_frameId-1}];
                                scatter3(xyz(moveID,1),xyz(moveID,2),xyz(moveID,3),45,(xyz(moveID,4)),'filled');
                            end
                            colormap(gca,"jet");
                            c = colorbar;
                            c.Label.String = 'velocity (m/s)';                        
                            grid on;
                            
                            xlim([-10 10])
                            ylim([1 maxRangeShow/2])
                            %zlim([-4 4])
                            zlim([-5 5])
                            xlabel('X (m)')
                            ylabel('y (m)')
                            zlabel('Z (m)')                        
                            
                            view(2)                        
                            title(' 3D point cloud');
                            
                            % plot range and azimuth heatmap
                            subplot(2,3,3)
                            % STATIC_ONLY: 1 = plot heatmap for zero-Doppler; 0 = plot heatmap for nonzero-Doppler
                            PLOT_ON_RA = 1;
                            STATIC_ONLY = 1; 
                            minRangeBinKeep = 1;
                            rightRangeBinDiscard = 1;
                            [mag_data_static(:,:,valid_obj_frameId), ...
                            mag_data_dynamic(:,:,valid_obj_frameId), ...
                            y_axis, x_axis] = plot_range_azimuth_2D(CFAR_DetectionObj.rangeBinSize, ...
                                                                    DopplerFFTOut,...
                                                                    length(genCalibrationMatrixObj.TxForMIMOProcess), ...
                                                                    length(genCalibrationMatrixObj.RxForMIMOProcess), ...
                                                                    CFAR_DetectionObj.antenna_azimuthonly, ...
                                                                    LOG_ON, ...
                                                                    STATIC_ONLY, ...
                                                                    PLOT_ON_RA, ...
                                                                    minRangeBinKeep, ...
                                                                    rightRangeBinDiscard);
                            title('range/azimuth heat map static objects')
                                
                            if PLOT_ON
                                subplot(2,3,4);
                                surf(y_axis, x_axis, (mag_data_static(:,:,valid_obj_frameId)).^0.4,'EdgeColor','none');
                                colormap(gca,"jet");
                                xlim([-10 10])
                                ylim([1 maxRangeShow/2])
                                view(2);
                                xlabel('meters');    
                                ylabel('meters');
                                title({'Static Range-Azimuth Heatmap',strcat('Current Total Frame Number = ', num2str(valid_obj_frameId))})
                                
                                subplot(2,3,5);
                                surf(y_axis, x_axis, (mag_data_dynamic(:,:,valid_obj_frameId)).^0.4,'EdgeColor','none');
                                colormap(gca,"jet");
                                xlim([-10 10])
                                ylim([1 maxRangeShow/2])
                                view(2);    
                                xlabel('meters');    
                                ylabel('meters');
                                title({'Dynamic HeatMap', strcat('curent Sweep angle= ',num2str(TxBF_Angle))});
                            end
                            
                            pause(0.1) 

                        end % PLOT_ON
            
                    end % if length(angleEst) > 0
                    
                    if ~mod(frameId,5)
                        % save Fig
                        temp = split(dataFolderName, '\');
                        angleName = temp(end-1);
                        temp(end)=[];
                        temp(end)=[];
                        tempStr = temp(1);
                        for iTemp = 2:length(temp)
                            tempStr = strcat(tempStr,'\', temp(iTemp));
                        end
                        tempStr = cell2mat(tempStr);
                        png_floder = strcat(tempStr, '_png\');
                        if ~exist(png_floder,'dir')
                            mkdir(png_floder)
                        end
                        p_file = strcat([png_floder, cell2mat(angleName), '_frame_',num2str(frameId), '.png']);
                        saveas(fig103, p_file, 'png');
                        %p_filefig = strcat([png_floder, cell2mat(angleName), '_frame_',num2str(frameId), '.fig']);
                        %saveas(fig103, p_filefig, 'fig');
                    end
                    
                    if LOG_ON
                        disp({strcat("Processing frameId >>>>> ", num2str(frameId)), ...
                            strcat("Processing valid_obj_frameId >>>>> ", num2str(valid_obj_frameId))});  
                    end
                end % PLOT_ON
                
                % return data
                gframe_obj{valid_obj_frameId}.rangeBins                     = gframe_obj_rangeBins(valid_obj_frameId, 1:length(obj_rangebin));
                gframe_obj{valid_obj_frameId}.estSNR                        = gframe_obj_estSNR(valid_obj_frameId, 1:length(obj_rangebin));
                gframe_obj{valid_obj_frameId}.estAngle                      = gframe_obj_estAngle(valid_obj_frameId, :, 1:length(obj_rangebin));                    
                gframe_obj{valid_obj_frameId}.RA_Spectrum                   = gframe_obj_RA_Spectrum(valid_obj_frameId,:);
                gframe_obj{valid_obj_frameId}.rangeResolution               = genCalibrationMatrixObj.rangeResolution;
                
                gframe_obj{valid_obj_frameId}.Range_Profile                 = rangeFFTOut(:,:,:);
                gframe_obj{valid_obj_frameId}.RangeDoppler_Profile          = DopplerFFTOut(:, [1:size(DopplerFFTOut,2)/2-1 size(DopplerFFTOut,2)/2+3:end], :);
                


                PLOT_ON_RA = 0;
                STATIC_ONLY = 0; 
                minRangeBinKeep = 1;
                rightRangeBinDiscard = 1;
                [mag_data_static(:,:,valid_obj_frameId), ...
                mag_data_dynamic(:,:,valid_obj_frameId), ...
                y_axis, x_axis] = plot_range_azimuth_2D(CFAR_DetectionObj.rangeBinSize, ...
                                                        DopplerFFTOut,...
                                                        length(genCalibrationMatrixObj.TxForMIMOProcess), ...
                                                        length(genCalibrationMatrixObj.RxForMIMOProcess), ...
                                                        CFAR_DetectionObj.antenna_azimuthonly, ...
                                                        LOG_ON, ...
                                                        STATIC_ONLY, ...
                                                        PLOT_ON_RA, ...
                                                        minRangeBinKeep, ...
                                                        rightRangeBinDiscard);
                

                gframe_obj{valid_obj_frameId}.RangeAngle_Dynamic_Profile    = mag_data_dynamic(:,:,valid_obj_frameId)';



            end % ~isempty(angleEst) 

        end % ~isempty(CFAR_detection_results)

    end % frameId


    %toc



end

