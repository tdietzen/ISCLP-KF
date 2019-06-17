function plotSpec(spec, magPowLin, xTickProp, yTickProp, cRange, plotColorbar)
% plotSpec(spec, magPowLin, xTickProp, yTickProp, cRange, plotColorbar)
% plots spectrogram.
%
% IN:
% spec              spectrogram data - freqbins x frames
% magPowLin         {'mag', 'pow', 'lin'}, scaling option, plot logarithmic magnitudes, logarithmic powers, or linear
% xTickProp         [xTick start, xTick resolution, xTick distance]
% xTickProp         [yTick start, yTick resolution, yTick distance]
% cRange            range of colorbar
% plotColorbar      plot colorbar if 1

% plot spectrogram
switch magPowLin
    case 'mag'
        imagesc(mag2db(abs(spec)), cRange); 
    case 'pow'
        imagesc(pow2db(abs(spec)), cRange); 
    case 'lin'
        imagesc(spec,              cRange);
    otherwise
        error('undefined scaling option.')   
end
set(gca,'Ydir','normal'); 
set(gca,'TickLength',[0 0]);

% add colorbar
if plotColorbar
    colorbar;
end

% add ticks
if ~isempty(xTickProp)
    xStartValue = xTickProp(1);
    xRes        = xTickProp(2);
    xTickDist   = xTickProp(3);  
    xTicks      = 1:xTickDist:size(spec,2);
    xTickLabels = cellstr(num2str((xStartValue + (xTicks-1)*xRes).'));
    set(gca,'xTick',xTicks)
    set(gca,'xTickLabel',xTickLabels);
else
    set(gca,'xTick',[])
end

if ~isempty(yTickProp)
    yStartValue = yTickProp(1);
    yRes        = yTickProp(2);
    ytickDist   = yTickProp(3);
    yTicks      = 1:ytickDist:size(spec,1);
    yTickLabels = cellstr(num2str((yStartValue + (yTicks-1)*yRes).'));
    set(gca,'yTick',yTicks);
    set(gca,'yTickLabel',yTickLabels);
else
    set(gca,'yTick',[]);
end

end