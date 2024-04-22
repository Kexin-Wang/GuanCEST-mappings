function [FitResult, FitParam] = PLOF(Offset, Saturation, FitParam)
% Fit the Z-spectrum using polynomial and Lorentzian line-shape fitting (PLOF) method

% If you have used this function in a scientific publication,we would
% appreciate citations to the following papers:

% [1]	Chen L, Wei Z, Cai S, Li Y, Liu G, Lu H, Weiss RG, van Zijl PCM, Xu
% J. High-resolution creatine mapping of mouse brain at 11.7 T using
% non-steady-state chemical exchange saturation transfer. NMR Biomed
% 2019:e4168. 
% [2]	Chen L, Barker PB, Weiss RG, van Zijl PCM, Xu J.
% Creatine and phosphocreatine mapping of mouse skeletal muscle by a
% polynomial and Lorentzian line-shape fitting CEST method. Magn Reson Med
% 2019;81(1):69-78. 
% [3]	Chen L, Zeng H, Xu X, Yadav NN, Cai S, Puts NA,
% Barker PB, Li T, Weiss RG, van Zijl PCM, Xu J. Investigation of the
% contribution of total creatine to the CEST Z-spectrum of brain using a
% knockout mouse model. NMR Biomed 2017;30(12):e3834.
% [4]   Wang K, Park S, Kamson DO, Li Y, Liu G, Xu J. 
%  Guanidinium and amide CEST mapping of human brain by high spectral resolution 
%  CEST at 3 T. Magn Reson Med. 2023 Jan;89(1):177-191. 

% Please contact Kexin Wang at kwang101@jhu.edu if you have any questions about the code. 

% UPDATED 11/14/2022:
%           update Rback in CurveFunction
% UPDATED 01/26/2023:
%           update x0, lb, ub for different B1
% UPDATED 11/17/2023:
%           update x0, lb, ub for B1=0.8uT for GRASE PD_CEST study

warning off
if size(Offset,1) == 1
    Offset = Offset';
end

if size(Saturation, 1) == 1
    Saturation = Saturation';
end

if (Saturation(1)) > 10
    error('Z-spectrum have to be 0-1, can not be percentage')
end


[Offset, Saturation] = Extract_Zspectrum(Offset, Saturation, FitParam.WholeRange);
[Offset_background, Saturation_background] = CutOff_Zspectrum(Offset, Saturation, FitParam.PeakRange);

FitResult.Offset = Offset;
FitResult.Saturation = Saturation;
FitResult.Offset_background = Offset_background;
FitResult.Saturation_background = Saturation_background;

%-----------------------------------------------%
%-------------- Two-step fitting ---------------%
%-----------------------------------------------%

FitResult.xindex = linspace(min(Offset), max(Offset), 100)';
%   First scale to hit 0, then scale to 3.5 ppm
if FitParam.satpwr == 0.3
        x0 = [3.2, .54, .87, -6.9, ...
                0.04, 1.50, 3.5, ...
                0.05, 1.31, 2.];
        lb = [1., 0.2, 0.1, -10,...
                7e-3, 0.4, 3.4, ...
                7e-3, 0.4, 1.9];
        ub = [5, 2, 2, 0., ...
                    1, 4, 3.7, ...
                    1, 4, 2.1];
elseif FitParam.satpwr == 0.4
        x0 = [10.34, 0.2945, .85, -7.8, ...
                0.0208, 0.64, 3.5, ...
                0.063, 2.14, 2.];
        lb = [1., 0.1, 0.1, -20,...
                7e-3, 0.4, 3.4, ...
                7e-3, 0.4, 1.9];
        ub = [15, 2, 2, 0., ...
                    1, 4, 3.6, ...
                    1, 4, 2.1];
elseif FitParam.satpwr == 0.6
        x0 = [9, 0.8, 0.8, -14, ...
                0.0356, 1.1093, 3.5, ...
                0.0851, 2.0328, 2.];
        lb = [1., .2, 0.5, -100,...
                7e-3, 0.4, 3.4, ...
                7e-3, 0.4, 1.9];
        ub = [200, 3.5, 1.8, 0., ...
                    .05, 4, 3.6, ...
                    .05, 4, 2.2];
% elseif FitParam.satpwr == 0.8
%         x0 = [3.7104, 1.0072, 1.0599, -14.1058, ...
%                 0.0346, 1.4423, 3.3242, ...
%                 0.0049, 0.6958, 2.1];
%         lb = [0, 0, 0, -100,...
%                 9e-3, 0.2, 3.3, ...
%                 9e-3, 0.2, 1.8];
%         ub = [10, 5, 5, -5, ...
%                     .08, 4, 3.7, ...
%                     .05, 4, 2.2];

% %   Fine-tune bg fitting and increase the lb for peak height and width of
% %   Guan
% elseif FitParam.satpwr == 0.8
%         x0 = [3.6, 0.96, 0.98, -12.3271, ...
%                 0.03, 0.66, 3.5, ...
%                 0.03, 0.1, 2.1];
%         lb = [0, 0, 0, -100,...
%                 .01, 0.04, 3.3, ...
%                 .01, 0.4, 1.8];
%         ub = [10, 5, 5, -5, ...
%                     .1, 4, 3.7, ...
%                     .1, 4, 2.2];

elseif FitParam.satpwr == 0.8
        x0 = [3, 0.60, 0.91, -14.07, ...
                0.03, 1.46, 3.35, ...
                0.07, 0.467, 2.1];
        lb = [0, 0, 0, -100,...
                9e-3, 0.02, 3.3, ...
                9e-3, 0.02, 1.8];
        ub = [1000, 5, 5, -5, ...
                    .08, 4, 3.7, ...
                    .05, 4, 2.2];

% %   used
% elseif FitParam.satpwr == 0.8
%         x0 = [3.1, 1.4, 2.6004, -17.3271, ...
%                 0.0799, 1.3894, 3.35, ...
%                 0.05, 1.27, 2.1];
%         lb = [0, 0, 0, -100,...
%                 9e-4, 0.02, 3.3, ...
%                 9e-4, 0.02, 1.8];
%         ub = [2000, 5, 5, -5, ...
%                     .08, 4, 3.7, ...
%                     .08, 4, 2.2];

        %   Previous for CrMap project
% elseif FitParam.satpwr == 0.8
%         x0 = [3.7104, 1.0072, 1.0599, -14.1058, ...
%                 0.0346, 1.4423, 3.3242, ...
%                 0.0049, 0.6958, 2.1];
%         lb = [0, 0, 0, -100,...
%                 9e-4, 0.02, 3.3, ...
%                 9e-4, 0.02, 1.8];
%         ub = [100, 5, 5, -5, ...
%                     .08, 4, 3.7, ...
%                     .05, 4, 2.2];
elseif FitParam.satpwr == 1.0
        x0 = [12.65, 0.7132, 1.2256, -19.9167, ...
                0.0607, 1.2324, 3.4216, ...
                0.1033, 1.6976, 2.];
        lb = [1., .2, 0.2, -100,...
                7e-3, 0.4, 3.4, ...
                7e-3, 0.4, 1.9];
        ub = [50, 5, 5, 0., ...
                    1, 4, 3.6, ...
                    1, 4, 2.1];
elseif FitParam.satpwr == 1.2
        x0 = [7.4, 1.2, 1.3, -30, ...
                0.048, 1.08, 3.4, ...
                0.057, 1.14, 2.];
        lb = [1., .2, 0.1, -100,...
                7e-3, 0.4, 3.4, ...
                7e-3, 0.4, 1.9];
        ub = [50, 5, 5, -10, ...
                    1, 4, 3.6, ...
                    1, 4, 2.1];
elseif FitParam.satpwr == 1.4
        x0 = [8.4, 1.3, 1.4, -32, ...
                0.05, 1.0, 3.5, ...
                0.06, 1.27, 2.];
        lb = [0., .2, 1e-4, -100,...
                7e-3, 0.4, 3.4, ...
                7e-3, 0.4, 1.9];
        ub = [50, 5, 5, -10, ...
                    1, 4, 3.6, ...
                    1, 4, 2.1];
elseif FitParam.satpwr == 1.6
        x0 = [9.1, 1.4, 1.5, -37, ...
                0.05, 1.07, 3.5, ...
                0.07, 1.3, 2.];
        lb = [.1, 1., 1, -1000,...
                7e-3, 0.4, 3.4, ...
                7e-3, 0.4, 1.9];
        ub = [50, 100, 100, -10, ...
                    1, 4, 3.6, ...
                    1, 4, 2.1];
end
%-------------background fitting --------------%

FitParam.peak = 0;
options = optimset('MaxFunEvals', 1e8,'TolFun', 1e-8,'TolX', 1e-8, 'Display',  'off' );
[P,SSE]=lsqcurvefit(@(p,xdata) CurveFunction(p,xdata,FitParam),x0, Offset_background, Saturation_background,lb,ub,options);
FitResult.Coefficents = P;

% fitting result of Z-spectrum and R-spectrum
[FitResult.Background, Rbak] = CurveFunction(FitResult.Coefficents, FitResult.xindex, FitParam);
Backgroundtmp = CurveFunction(FitResult.Coefficents, Offset, FitParam);
FitResult.DeltaZspectrum = Backgroundtmp - Saturation; 

% R^2 for the evaluation of fitting
FitResult.RsquareBG = 1 - goodnessOfFit(CurveFunction(FitResult.Coefficents, Offset_background, FitParam), Saturation_background, 'NMSE');

%   fitted background at 3.5 ppm
tmplist = abs(FitResult.xindex - 3.5);
[res,  idx] = min(tmplist);
FitResult.MTbackground = FitResult.Background(idx); 

%-------------CEST all-peak fitting --------------%
FitParam.peak = 1;

% fix the background parameters
lb(1: 4) = FitResult.Coefficents(1 : 4); 
ub(1: 4) = FitResult.Coefficents(1 : 4); 

[FitResult.Coefficents, resnorm] = lsqcurvefit(@CurveFunction, x0, Offset, Saturation, lb, ub, options, FitParam);
[FitResult.Curve, FitResult.Rcurve] = CurveFunction(FitResult.Coefficents, FitResult.xindex, FitParam);

FitResult.DeltaFitZ = FitResult.Background-FitResult.Curve;
FitResult.DeltaFitR = FitResult.Rcurve - Rbak;


%------------ Calculate deltaZ and assign the Coefficents ---------------%

FitResult.Rpeak1 = FitResult.Coefficents(5);
FitResult.FitPeakOffset1 = FitResult.Coefficents(7);
FitResult.Rpeak2 = FitResult.Coefficents(8);
FitResult.FitPeakOffset2 = FitResult.Coefficents(10); 

FitResult.MT = FitResult.Coefficents(1); 
FitResult.RsquareAll = 1 - goodnessOfFit(CurveFunction(FitResult.Coefficents, Offset, FitParam), Saturation, 'NMSE');


%------------ Calculate DeltaZ-spectrum ---------------%

% amide peak
FitParam.peak = 2;

[Ztmp Rtmp]= CurveFunction(FitResult.Coefficents, FitResult.xindex, FitParam); 
FitResult.ZamideFit =FitResult.Background -Ztmp;
[FitResult.DeltaZpeak1,  idx]= max(FitResult.ZamideFit);
FitResult.DeltaZpeak1Offset = FitResult.xindex(idx); 
FitResult.RamideFit =Rtmp - Rbak;
[FitResult.DeltaRpeak1, idx] = max(FitResult.RamideFit);
FitResult.DeltaRpeak1Offset = FitResult.xindex(idx); 

% guan peak
FitParam.peak = 3;

[Ztmp Rtmp] = CurveFunction(FitResult.Coefficents, FitResult.xindex, FitParam); 
FitResult.ZguanFit =FitResult.Background -Ztmp;
[FitResult.DeltaZpeak2,  idx] = max(FitResult.ZguanFit);
FitResult.DeltaZpeak2Offset = FitResult.xindex(idx); 
FitResult.RguanFit =Rtmp - Rbak;
[FitResult.DeltaRpeak2,  idx] = max(FitResult.RguanFit);
FitResult.DeltaRpeak2Offset = FitResult.xindex(idx); 

%------------ modify FitResult.SNR ---------------%
if FitParam.CalSNR == 1
    sample = FitParam.SNRimg;
    SNRmask = FitParam.SNRmask;

    signal = mean(sample(SNRmask));

    [M, N] = size(sample);
    maskCorners = zeros(M, N);
    t = 16; % ensure the sum of 4 squares in the corner >= 1000 pixels
    maskCorners(1:t, 1:t) = 1;
    maskCorners(1:t, N - t + 1: N) = 1;
    maskCorners(M - t + 1: M, 1:t) = 1;
    maskCorners(M - t + 1: M, N - t + 1: N) = 1;
    noise = std(sample(logical(maskCorners)));
    FitResult.SNR = signal/noise;
end


if FitParam.ifshowimage == 1
    PlotFitResult(FitResult,FitParam);
end
end