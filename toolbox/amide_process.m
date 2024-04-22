function [displayimg_origin, FitResult] = amide_process(mask,cestimgs,Offsets,T1_map, FitParam)
%   INPUT:
%           mask: ROI of the brain
%           cestimgs: CEST images
%           Offsets: frequency list
%           T1_map: T1 map
%           stp: saturation power [uT]
%           FitParam: fitting parameters
%   OUTPUT:
%           displayimg_origin: CEST images
%           cestinspect: Z spectrum
%           FitResult: save all the GuanCEST and amideCEST maps, etc.

%   UPDATED 2022/10/19:
%           Change the fitting peak range to (0.9, 5.0) from (1.5, 5.0) ppm
%   UPDATED 2023/02/24:
%           1. Add T1 map check 
%           2. Add the scaling of the Z-spectrum
%           3. Add CrMap
%   UPDATE 2023/03/08:
%           1. Add dZex and dZfit
%   UPDATED 2023/05/09:
%           1. Resize the cestimgs, if its size mismatched that of the mask
%   UPDATED 2023/08/09:
%           1. Add the replacement of NaN and inf in the cest images before
%           MLSVD
%   UPDATED 2023/08/11:
%           1. Add the conditional statement to avoid the case that cest
%           images are resliced by SPM and present smaller regions than the
%           masks
%   UPDATED 2023/11/29:
%           1. Remove the denoise, and generally keep all the pixels.
%           2. Increase the mlsvd core size to have a rougher smoothing and
%           redunce the number of bad points.
%   UPDATED 2023/12/18:
%           1. Use previous MLSVD core size: [Ue, Se, ~] = mlsvd(cestimgs,[48 48 10]);

    if size(mask, 1) ~= size(cestimgs, 1)
            wanted_size = [size(mask, 1), size(mask, 2), size(cestimgs, 3)];
            cestimgs = imresizen(cestimgs,wanted_size./size(cestimgs));
    end

    IMGtmp_origin = cestimgs;
    displayimg_origin = IMGtmp_origin(:,:,1); 
    displayimg_origin(~mask) = 0;
    scaler = mean(displayimg_origin(mask));
    displayimg_origin = displayimg_origin./scaler;
    %   For the scaling of the Z-spectrum
%     xq = flipud((-1.:0.01:8.)');
%     idx_list=1:10:length(xq);
%     idx_list=idx_list';
%     fullppm=xq(idx_list);

    % compute the MLSVD of the tensor
    %   Replace NaN and inf with 0 before MLSVD
    cestimgs(~isfinite(cestimgs)) = 0;
    %   All the NaN/Inf/truncated regions by SPM will be removed from the
    %   mask and the following analysis
    mask = logical(mask .* logical(squeeze(cestimgs(:,:,2))));
    % [Ue, Se, ~] = mlsvd(cestimgs,[72 72 18]);
    [Ue, Se, ~] = mlsvd(cestimgs,[48 48 10]);
    IMGtmp_process = lmlragen(Ue, Se);


    %refine MASK image
    displayimg_process = IMGtmp_process(:,:,1);
    displayimg_process(~mask) = 0;
    NXALL = size(IMGtmp_process,1);
    NYALL = size(IMGtmp_process,2);
    % for idx = 1:NXALL
    %    for idy = 1:NYALL
    %         if (mask(idx,idy))
    %             % if (displayimg_process(idx,idy) < 0.7*mean2(displayimg_process(mask)))
    %                 mask(idx,idy) = false;
    %             end
    %         end
    %    end
    % end


     CESTimg = IMGtmp_process(:,:,:)./IMGtmp_process(:,:,2);

     %remove first two M0 points
     FreqPPM = Offsets(3:end);
     CESTimg = CESTimg(:,:,3:end);
     
%      %get whole maks CEST Z-spectrum for inspection
%      for i = 1:size(FreqPPM)
%          cestinspect(:,1) = FreqPPM;
%          tmp = CESTimg(:,:,i);
%          cestinspect(i,2) = mean(tmp(mask));
%      end
     
     % fit each point
     idxall=0;
    
     Zamide = zeros(NXALL, NYALL);
     Zguan = zeros(NXALL, NYALL);
     Ramide = zeros(NXALL, NYALL);
     Rguan = zeros(NXALL, NYALL);
     MTbg = zeros(NXALL, NYALL);
     NOE = zeros(NXALL, NYALL);
%      Zcr = zeros(NXALL, NYALL);
     T1flag = mean(T1_map(mask));
     dZfit = [];
     dZex = [];
     Zfit = [];
     Zex = [];
     flist = [];
     snr = [];
     Zbg = [];
     for idx=1:NXALL
        for idy=1:NYALL
            idxall=idxall+1;

            if (mask(idx,idy))
                          
                Z_spectrum = CESTimg(idx,idy,:);
                Z_spectrum = squeeze(Z_spectrum);
        
                % T1 correction
                %   T1 check, if T1 is in the unit of ms, change it to s
                if T1_map(idx,idy) == 0
                    FitParam.R1 = 1;
                elseif T1flag > 10
                    FitParam.R1 = 1/(T1_map(idx,idy)/1e3);
                else
                    FitParam.R1 = 1/T1_map(idx,idy);
                end 
                                %   Scaled to [0, 1] by (Z - Zmin)./(1 - Zmin)
                    [Zmin, min_idx] = min(Z_spectrum);
                    z_tmp = (Z_spectrum - Zmin)./(1 - Zmin);
                % B0 shift
                FitParam.WholeRange = [-2, 2];
                B0Shift = WASSR(FreqPPM, z_tmp, FitParam);
                
                %%   Scale the Z-spectrum to hit zero and have a standard
%                 %   amide-Z
%                 %   Interpolation and scaled to 0
%                 z_tmp = interp1(FreqPPM - B0Shift, Z_spectrum,xq,'spline');
                
%                 %   For debug
%                 if idy == 44    %   idx ==32 
%                     FitParam.ifshowimage = 1; 
%                     figure();
%                     plot(FreqPPM, Z_spectrum);
% %                     hold on;
% %                     plot(xq, z_tmp);
% %                     hold off;
%                 else
%                     FitParam.ifshowimage = 0; 
%                 end



                    %   NOE
                    [tmp, idx1] = min(abs(FreqPPM + 10));
                    [tmp, idx2] = min(abs(FreqPPM + 3.5));
                    NOE(idx,idy) =100*(z_tmp(idx1) - z_tmp(idx2));
                
%                 %   Scaled to [0, 1] by (Z - Zmin)./(1 - Zmin)
%                 [Zmin, min_idx] = min(z_tmp);
%                 z_tmp = (z_tmp - Zmin)./(1 - Zmin);
%                 z_tmp = z_tmp(idx_list);
                %   Scaled for the amide
%                 amideValue = z_tmp(amide_idx);
%                 z_tmp = z_tmp./amideValue.*scale.standard_amideValue;
                                
                %% Fit PLOT 
                FitParam.WholeRange = [0, 8];    % CEST peak parameters
                if z_tmp(1) > 3
                    mask(idx,idy) = 0;
                    continue
                end
                [FitResult,FitParam] = PLOF(FreqPPM - B0Shift, z_tmp, FitParam);
                % if size(FitResult.DeltaZspectrum, 1) > 40
                        Zamide(idx,idy) =100*FitResult.DeltaZpeak1;
                        Ramide(idx,idy) = 1000*FitResult.DeltaRpeak1;
                        Zguan(idx,idy) =100*FitResult.DeltaZpeak2;
                        Rguan(idx,idy) = 1000*FitResult.DeltaRpeak2;
                        MTbg(idx,idy) = 100 * (1 - FitResult.MTbackground);
        %                 Zcr(idx,idy) = Zguan(idx,idy) - scale.aar * Zamide(idx,idy);
                        %   (1:75) may be changed if different crlist is used
                        % dZex = [dZex, FitResult.DeltaZspectrum(1:40)];
                        dZfit = [dZfit, FitResult.DeltaFitZ];
                        Zfit = [Zfit, FitResult.Curve];
                        snr = [snr, FitResult.SNR];
                        % Zex = [Zex, FitResult.Saturation(1:40)];
                        % flist = [flist, FitResult.Offset(1:40)];      
                        Zbg = [Zbg, FitResult.Background];
                % else
                %         mask(idx,idy) = false;
                % end

                                               
            end
            waitbar(idxall/(NXALL*NYALL));
        end
     end
    % xindex = FitResult.xindex;
    % dZ_Offset = FitResult.Offset;
    clear FitResult;
    % FitResult.xindex = xindex;
    % FitResult.Offset = dZ_Offset;
    FitResult.ZamideMap = Zamide;
    FitResult.ZguanMap = Zguan;
    FitResult.RamideMap = Ramide;
    FitResult.RguanMap = Rguan;
    FitResult.MTbgMap = MTbg;
    FitResult.rNOEMap = NOE;
%     FitResult.ZcrMap = Zcr;
    % FitResult.dZex = mean(dZex, 2);
    FitResult.dZfit = mean(dZfit, 2);
    FitResult.Zfit = mean(Zfit, 2);
    FitResult.mask = mask;
    FitResult.SNR = mean(snr);
    % FitResult.Zex = mean(Zex, 2);
    % FitResult.flist = mean(flist, 2);
    FitResult.Zbg = mean(Zbg, 2);
end