function [] = drawShadedRect(rectWidth, rectHeight,  xCtr , yCtr, angle, lineColor, isOutlier, ax)

% set up outline of unrotated rectangle
xStart = xCtr - (rectWidth-1)/2;
yStart = yCtr - (rectHeight-1)/2;
xEnd = xCtr + (rectWidth-1)/2;
yEnd = yCtr + (rectHeight-1)/2;
xv = [xStart xEnd xEnd xStart xStart  xCtr xCtr ];% xEnd];
yv = [yStart yStart yEnd yEnd yStart yStart yCtr];% yCtr];


% rotate the rectangle around its center
XYOrig(1,:)=xv-xCtr;XYOrig(2,:)=yv-yCtr;
rotXYOrig=[cosd(angle) -sind(angle);sind(angle) cosd(angle)]*XYOrig;
rotXY(1,:)=rotXYOrig(1,:)+xCtr;
rotXY(2,:)=rotXYOrig(2,:)+yCtr;



if 1
    % % plot a colored patch
    patch(rotXY(1,1:end-2),rotXY(2,1:end-2),lineColor,'parent',ax,'facealpha',.25,'edgecolor','none');
    % plot outline around patch
end


set(ax,'ydir','reverse')


hold(ax, 'on');
if ~isOutlier
    plot(ax,rotXY(1,end-1:end),rotXY(2,end-1:end),'color',lineColor,'linewidth',1);
else
    plot(ax,rotXY(1,end-1:end),rotXY(2,end-1:end),'-','color',lineColor,'linewidth',10);
end

