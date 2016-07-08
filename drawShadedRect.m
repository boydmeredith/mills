function [] = drawShadedRect(rectWidth, rectHeight,  xCtr , yCtr, angle, patchColor,isOutlier, ax)

% set up outline of unrotated rectangle
xStart = xCtr - (rectWidth-1)/2;
yStart = yCtr - (rectHeight-1)/2;
xEnd = xCtr + (rectWidth-1)/2;
yEnd = yCtr + (rectHeight-1)/2;
xv = [xStart xEnd xEnd xStart xStart  xCtr xCtr xEnd];
yv = [yStart yStart yEnd yEnd yStart yStart yCtr yCtr];


% rotate the rectangle around its center
XYOrig(1,:)=xv-xCtr;XYOrig(2,:)=yv-yCtr;
rotXYOrig=[cosd(angle) -sind(angle);sind(angle) cosd(angle)]*XYOrig;
rotXY(1,:)=rotXYOrig(1,:)+xCtr;
rotXY(2,:)=rotXYOrig(2,:)+yCtr;

% plot a colored patch
patch(rotXY(1,1:end-3),rotXY(2,1:end-3),patchColor,'parent',ax);
% plot outline around patch
hold(ax, 'on');
if ~isOutlier
    plot(ax,rotXY(1,:),rotXY(2,:),'k');
else
    plot(ax,rotXY(1,:),rotXY(2,:),'c');
end
set(ax,'ydir','reverse')
colormap(ax, colormapRedBlue)
colorbar(ax);
