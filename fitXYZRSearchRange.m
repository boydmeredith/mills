function [xyzrSearchRange, outliersXY] = fitXYZRSearchRange(x,y,z,r,pRes,doMakePlot)
% function [xyzrSearchRange, outliersXY] = fitXYZRSearchRange(x,y,z,r,pRes)
%
% estimate fits for x, y, z, r. use x, y to find outliers, since they
% should fit nicely in a grid.
%
% input:
% x, y, z, r
% pRes
%
% output: 
% xyzrSearchRange
% outliersXY

if nargin < 6, doMakePlot = false, end;

xyzrSearchRange = nan(4,pRes.nBlockSpan^2);
% create a grid of block i,j indices
[xx, yy] = meshgrid(1:pRes.nBlockSpan,1:pRes.nBlockSpan);
% robustly fit xs and ys to plane
[fX, ~, outX]  = fit([xx(:), yy(:)], reshape(x, pRes.nBlockSpan^2, 1), 'poly11','Robust','Bisquare');
[fY, ~, outY]  = fit([xx(:), yy(:)], reshape(y, pRes.nBlockSpan^2, 1), 'poly11','Robust','Bisquare');
% find x and y outliers
outliersX = abs(outX.residuals) > pRes.nRSTD*robustSTD(outX.residuals);
outliersY = abs(outY.residuals) > pRes.nRSTD*robustSTD(outY.residuals);
outliersXY = outliersX' | outliersY';
% put the xy fit values into the search range
xyzrSearchRange(1,:) = fX(xx(:),yy(:));
xyzrSearchRange(2,:) = fY(xx(:),yy(:));
% if there are zs to fit, fit them, excluding the XY outliers
if ~isempty(z)
    z = reshape(z, pRes.nBlockSpan^2, 1);
    [fZ, ~, ~] = fit([xyzrSearchRange(1,~outliersXY)', xyzrSearchRange(2,~outliersXY)'], z(~outliersXY), 'loess','Robust','off');
    % round the fit values to the nearest integer so we can use them as neighborhood centers
    xyzrSearchRange(3,:) = round(fZ(xyzrSearchRange(1,:),xyzrSearchRange(2,:)));
end
%  plot(fZ,[xx(:),yy(:)],fZ(xx(:),yy(:)))
% if there are rs to fit, fit them, excluding the XY outliers
if ~isempty(r)
    r = reshape(r, pRes.nBlockSpan^2, 1);
    [fR, ~, outR] = fit([xyzrSearchRange(1,~outliersXY)',xyzrSearchRange(2,~outliersXY)'], r(~outliersXY), 'poly11','Robust','off');
    % make sure that the rotation angles are divisible by the fine
    % rotation angle step size parameter
    rFits = fR(xyzrSearchRange(1,:)',xyzrSearchRange(2,:)'); 
    xyzrSearchRange(4,:) = round(rFits/pRes.fineRotStepSz) * pRes.fineRotStepSz);
end
%  plot(fR,[xx(:),yy(:)],fR(xx(:),yy(:)))


if ~isempty(pRes.searchRangeFigName) & doMakePlot
    searchRangeFig = figure('visible', pRes.showFigs);
    
    [~, pfM]      = makeSubplots(searchRangeFig, 2, 5, .1, .5, [0 0 1 1]);
    colormap(searchRangeFig, colormapRedBlue);
    
    imagesc([reshape(x,pRes.nBlockSpan,pRes.nBlockSpan)...
        reshape(xyzrSearchRange(1,:),pRes.nBlockSpan,pRes.nBlockSpan)],'parent',pfM(1,1));
    title(pfM(1,1),'X');     colorbar(pfM(1,1))
    axis(pfM(1,1),'image');
    
    imagesc(reshape(outliersX,pRes.nBlockSpan,pRes.nBlockSpan),'parent',pfM(2,1));
    axis(pfM(2,1),'image');
    title(pfM(2,1),'X outliers')
    
    hist(pfM(3,1), outX.residuals,1000);
    hold(pfM(3,1), 'on')
    plot(pfM(3,1), pRes.nRSTD*[robustSTD(outX.residuals) robustSTD(outX.residuals)],get(pfM(3,1),'ylim'),'r',...
        -pRes.nRSTD*[robustSTD(outX.residuals) robustSTD(outX.residuals)],get(pfM(3,1),'ylim'),'r')
    set(pfM(3,1),'ycolor','w','tickdir','out');box(pfM(3,1),'off');title(pfM(3,1),'X fit residuals')
    
    imagesc([reshape(y,pRes.nBlockSpan,pRes.nBlockSpan)...
        reshape(xyzrSearchRange(2,:),pRes.nBlockSpan,pRes.nBlockSpan)],'parent',pfM(1,2));
    title(pfM(1,2),'Y'); colorbar(pfM(1,2))
    axis(pfM(1,2),'image');
    
    imagesc(reshape(outliersY,pRes.nBlockSpan,pRes.nBlockSpan),'parent',pfM(2,2));
    axis(pfM(2,2),'image');
    title(pfM(2,1),'Y outliers')
    
    
    hist(pfM(3,2), outY.residuals,1000);
    hold(pfM(3,2), 'on')
    plot(pfM(3,2), pRes.nRSTD*[robustSTD(outY.residuals) robustSTD(outY.residuals)],get(pfM(3,2),'ylim'),'r',...
        -pRes.nRSTD*[robustSTD(outY.residuals) robustSTD(outY.residuals)],get(pfM(3,2),'ylim'),'r')
    set(pfM(3,2),'ycolor','w','tickdir','out');box(pfM(3,2),'off');title(pfM(3,2),'Y fit residuals')
    
    
    imagesc([reshape(z,pRes.nBlockSpan,pRes.nBlockSpan)...
        reshape(xyzrSearchRange(3,:),pRes.nBlockSpan,pRes.nBlockSpan)],'parent',pfM(4,1));
    title(pfM(4,1),'Z'); colorbar(pfM(4,1))
    axis(pfM(4,1),'image');
    
    imagesc([reshape(r,pRes.nBlockSpan,pRes.nBlockSpan)...
        reshape(xyzrSearchRange(4,:),pRes.nBlockSpan,pRes.nBlockSpan)],'parent',pfM(4,2));
    title(pfM(4,2),'R'); colorbar(pfM(4,2))
    axis(pfM(4,2),'image');
    
    plot(fZ,[xyzrSearchRange(1,:)',xyzrSearchRange(2,:)'],z,'parent',pfM(5,1))
    plot(fR,[xyzrSearchRange(1,:)',xyzrSearchRange(2,:)'],r,'parent',pfM(5,2))
    
    set(searchRangeFig, 'position',[50 50 1000 600], 'paperpositionmode','manual',...
        'paperunits','inches','paperposition',[0 0 8.5 11]);
    saveas(searchRangeFig, fullfile(pRes.corrDir, pRes.searchRangeFigName));
    %close(searchRangeFig)
end
end