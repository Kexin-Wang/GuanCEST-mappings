function [Z, R] = CurveFunction(x, xdata, FitParam)
%------------ Input------------%
%
%   x is the list of Lorentzian-peak related fitting parameters: 
%       x(1~4) for background;   
%       x(5~7) for peak_1;   x(8~10) for peak_2; ...
%
%   x =  [Amp_background, width_background, Coefficient_of_Zero_Order_background, Coefficient_of_First_Order_background, ...
%           Rexch_peak_1, width_peak_1, center_peak_1, ...
%           Rexch_peak_2, width_peak_2, center_peak_2, ...]
%
%   xdata is the list of offsets, y value (signal Z) is in the FitParam as
%   FitParam.Satruration
%
%   peak = FitParam.peak;
%   peak = 0: only fit background;
%           = 1: all peaks are used for the fitting;
%           = 2: only backgound and peak_1 are used for the fitting;
%           = 3: only backgound and peak_2 are used for the fitting;
%           ...
%
%------------ Output------------%
% 
%   Z is the Z-spectrum.


% Eq. (6): R_back = ; C_0/1/2/3 becomes x(4/5/6/7)
% peakoffset = FitParam.PeakOffset;
% Rbak = x(1)*x(2)^2./( x(2)^2 + 4*(xdata ).^2 ) + x(3) + x(4)*0.001*(xdata - peakoffset);
Rbak = x(1)*x(2)^2./( x(2)^2 + 4*(xdata ).^2 ) + x(3) + x(4)*0.001*xdata;

% Eq. (5): R = R_exch * (w/2) ^2 over (w/2) ^2 + (delta_w - delta_w_exch)^2
Ramide = x(5)*x(6)^2./( x(6)^2 + 4*(xdata - x(7)).^2 );
Rguan = x(8)*x(9)^2./( x(9)^2 + 4*(xdata - x(10)).^2 );

 % cos^2(theta)
satHz = FitParam.satpwr*42.58;
cos2thet = 1 - satHz^2./(satHz^2+ (FitParam.Magfield*xdata).^2);

peak = FitParam.peak;
if peak == 0 
    Rall = Rbak;
elseif peak == 1 
    Rall = Rbak + Ramide + Rguan;
elseif peak == 2
    Rall = Rbak + Ramide;
elseif peak == 3
    Rall = Rbak + Rguan;
end

% Output of Z-spectrum 
Z = (1 - cos2thet*FitParam.R1./Rall).* exp(-Rall * FitParam.tsat) + cos2thet * FitParam.R1./Rall;
R = Rall;
end