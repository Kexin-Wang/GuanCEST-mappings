function B0 = WASSR(Offset, Saturation, FitParam)

warning off
if size(Offset, 1) == 1
    Offset = Offset';
end

if size(Saturation, 1) == 1
    Saturation = Saturation';
end

% if (Saturation(1)) > 10
%     error('Z-spectrum have to be 0-1, can not be percentage')
% end

% B0 = 0;
[Offset, Saturation] = Extract_Zspectrum(Offset, Saturation, FitParam.WholeRange);
p = polyfit(Offset, Saturation, 6);
x1 = FitParam.WholeRange(1): 0.01: FitParam.WholeRange(2);
y1 = polyval(p, x1);
[tmp, idx] = min(y1);
B0 = x1(idx);

if FitParam.ifshowimage == 1
    figure;
    drawnow;
     hold on;
    plot(Offset, Saturation,'o');
    plot(x1, y1, '-');
   
   hold off;
end

end
