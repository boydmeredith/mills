function [xyzrSearchRange, outliersXY] = fitXYZRSearchRange(x,y,z,r,pRes,doMakePlot,doSavePlot)
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

if nargin < 6, doMakePlot = false; end;
if nargin < 7, doSavePlot = true; end;

xyzrSearchRange = nan(4,pRes.nBlockSpan^2);
% create a grid of block i,j indices
[xx, yy] = meshgrid(1:pRes.nBlockSpan,1:pRes.nBlockSpan);
% robustly fit xs and ys to plane
[fX, ~, outX]  = fit([xx(:), yy(:)], reshape(x, pRes.nBlockSpan^2, 1), 'poly11','Robust','Bisquare');
[fY, ~, outY]  = fit([xx(:), yy(:)], reshape(y, pRes.nBlockSpan^2, 1), 'poly11','Robust','Bisquare');
% find x and y outliers
xOutRadius = max(pRes.nRSTD*robustSTD(outX.residuals),pRes.xRadiusMin) ;
yOutRadius = max(pRes.nRSTD*robustSTD(outY.residuals),pRes.yRadiusMin) ;
outliersX = abs(outX.residuals) > xOutRadius;
outliersY = abs(outY.residuals) > yOutRadius;
outliersXY = outliersX' | outliersY';
% put the xy fit values into the search range
xyzrSearchRange(1,:) = fX(xx(:),yy(:));
xyzrSearchRange(2,:) = fY(xx(:),yy(:));
% if there are zs to fit, fit them, excluding the XY outliers
if ~isempty(z)
    z = reshape(z, pRes.nBlockSpan^2, 1);
    [fZ, ~, outZ] = fit([xyzrSearchRange(1,~outliersXY)', xyzrSearchRange(2,~outliersXY)'], z(~outliersXY), 'loess','Robust','off');
    % round the fit values to the nearest integer so we can use them as neighborhood centers
    zFits = fZ(xyzrSearchRange(1,:),xyzrSearchRange(2,:))';
    xyzrSearchRange(3,:) = round(zFits);
end
%  plot(fZ,[xx(:),yy(:)],fZ(xx(:),yy(:)))
% if there are rs to fit, fit them, excluding the XY outliers
if ~isempty(r)
    r = reshape(r, pRes.nBlockSpan^2, 1);
    [fR, ~, outR] = fit([xyzrSearchRange(1,~outliersXY)',xyzrSearchRange(2,~outliersXY)'], r(~outliersXY), 'poly11','Robust','off');
    % make sure that the rotation angles are divisible by the fine
    % rotation angle step size parameter
    rFits = fR(xyzrSearchRange(1,:)',xyzrSearchRange(2,:)'); 
    xyzrSearchRange(4,:) = round(rFits/pRes.fineRotStepSz) * pRes.fineRotStepSz;
end
%  plot(fR,[xx(:),yy(:)],fR(xx(:),yy(:)))


if ~isempty(pRes.searchRangeFigName) && doMakePlot
    searchRangeFig = figure('visible', pRes.showFigs);
    
    [flatP, flatM]      = makeSubplots(searchRangeFig, 2, 4, .2, .4, [0 0 .6 1]);
    [surfP, surfM]      = makeSubplots(searchRangeFig, 1, 2, .2, 0, [.6 0 .4 1]);
    colormap(searchRangeFig, colormapRedBlue);
    
    imagesc([reshape(x,pRes.nBlockSpan,pRes.nBlockSpan)...
        reshape(xyzrSearchRange(1,:),pRes.nBlockSpan,pRes.nBlockSpan)],'parent',flatM(1,1));
    title(flatM(1,1),'X');     colorbar('peer',flatM(1,1))
    axis(flatM(1,1),'image');
    
    imagesc(reshape(outliersX,pRes.nBlockSpan,pRes.nBlockSpan),'parent',flatM(2,1));
    axis(flatM(2,1),'image');
    title(flatM(2,1),'X outliers')
    
    hist(flatM(3,1), outX.residuals,1000);
    hold(flatM(3,1), 'on')
    plot(flatM(3,1), [xOutRadius xOutRadius],get(flatM(3,1),'ylim'),'r',...
        -[xOutRadius xOutRadius],get(flatM(3,1),'ylim'),'r')
    set(flatM(3,1),'ycolor','w','tickdir','out');box(flatM(3,1),'off');title(flatM(3,1),'X fit residuals')
    
    imagesc([reshape(y,pRes.nBlockSpan,pRes.nBlockSpan)...
        reshape(xyzrSearchRange(2,:),pRes.nBlockSpan,pRes.nBlockSpan)],'parent',flatM(1,2));
    title(flatM(1,2),'Y'); colorbar('peer',flatM(1,2))
    axis(flatM(1,2),'image');
    
    imagesc(reshape(outliersY,pRes.nBlockSpan,pRes.nBlockSpan),'parent',flatM(2,2));
    axis(flatM(2,2),'image');
    title(flatM(2,2),'Y outliers')
    title(flatM(1,2),'X outliers')

    
    hist(flatM(3,2), outY.residuals,1000);
    hold(flatM(3,2), 'on')
    plot(flatM(3,2), [yOutRadius yOutRadius],get(flatM(3,2),'ylim'),'r',...
        -[yOutRadius yOutRadius],get(flatM(3,2),'ylim'),'r')
    set(flatM(3,2),'ycolor','w','tickdir','out');box(flatM(3,2),'off');title(flatM(3,2),'Y fit residuals')
    
    
    imagesc([reshape(z,pRes.nBlockSpan,pRes.nBlockSpan)...
        reshape(xyzrSearchRange(3,:),pRes.nBlockSpan,pRes.nBlockSpan)],'parent',flatM(4,1));
    title(flatM(4,1),'Z'); colorbar('peer',flatM(4,1))
    axis(flatM(4,1),'image');
    
    imagesc([reshape(r,pRes.nBlockSpan,pRes.nBlockSpan)...
        reshape(xyzrSearchRange(4,:),pRes.nBlockSpan,pRes.nBlockSpan)],'parent',flatM(4,2));
    title(flatM(4,2),'R'); colorbar('peer',flatM(4,2))
    axis(flatM(4,2),'image');
    
    plot(fZ,[xyzrSearchRange(1,:)',xyzrSearchRange(2,:)'],z,'parent',surfM(1,1));
    plot(fZ,[xyzrSearchRange(1,:)',xyzrSearchRange(2,:)'],z,'parent',surfM(1,1));
    hold(surfM(1,1),'on');
    plot3(repmat(xyzrSearchRange(1,:),2,1), repmat(xyzrSearchRange(2,:),2,1), [zFits z]',...
        'k','linewidth',5,'parent',surfM(1,1))
    shading(surfM(1,1),'flat')
    axis(surfM(1,1),'square');
%     plot3(surfM(1,1), xyzrSearchRange(1,:),xyzrSearchRange(2,:),xyzrSearchRange(3,:)
    set(surfM(1,1),'view',[-40.0000   16.4000],'ydir','reverse','xlim',[0 512], 'ylim', [0 512])
    
    plot(fR,[xyzrSearchRange(1,:)',xyzrSearchRange(2,:)'],r,'parent',surfM(2,1));
    plot(fR,[xyzrSearchRange(1,:)',xyzrSearchRange(2,:)'],r,'parent',surfM(2,1));
    hold(surfM(2,1),'on')
    shading(surfM(2,1),'flat')
    plot3(repmat(xyzrSearchRange(1,:),2,1), repmat(xyzrSearchRange(2,:),2,1), [rFits r]',...
        'k','linewidth',5,'parent',surfM(2,1))
    xlabel(surfM(1,1),'x');
    ylabel(surfM(1,1),'y');
    zlabel(surfM(1,1),'z');
    axis(surfM(2,1),'square');
    set(surfM(2,1),'view',[-40.0000   16.4000],'ydir','reverse','xlim',[0 512], 'ylim', [0 512])
    xlabel(surfM(2,1),'x');
    ylabel(surfM(2,1),'y');
    zlabel(surfM(2,1),'rot angle');
    
    %%
%     plot(fZ,[xyzrSearchRange(1,:)',xyzrSearchRange(2,:)'],z,'style','residuals','parent',surfM(2,1));
%     axis(surfM(2,1),'square');
%     hold(surfM(2,1),'on');
% %     plot3(surfM(1,1), xyzrSearchRange(1,:),xyzrSearchRange(2,:),xyzrSearchRange(3,:)
%     set(surfM(2,1),'view',[-40.0000   16.4000],'ydir','reverse')
%     
%     plot(fR,[xyzrSearchRange(1,:)',xyzrSearchRange(2,:)'],r,'style','residuals','parent',surfM(2,2));
%     xlabel(surfM(2,1),'x');
%     ylabel(surfM(2,1),'y');
%     zlabel(surfM(2,1),'z');
%     axis(surfM(2,2),'square');
%  
%     set(surfM(2,2),'view',[-40.0000   16.4000],'ydir','reverse')
%     xlabel(surfM(2,2),'x');
%     ylabel(surfM(2,2),'y');
%     zlabel(surfM(2,2),'rot angle');
    
       set(searchRangeFig, 'position',[53 5 1220 700], 'paperpositionmode','manual',...
        'paperunits','inches','paperposition',[0 0 8.5 11],'paperorientation','landscape');
    
    if doSavePlot
    saveas(searchRangeFig, fullfile(pRes.corrDir, pRes.searchRangeFigName));
    end
    %close(searchRangeFig)
end
end
