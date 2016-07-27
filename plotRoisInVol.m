function  plotRoisInVol(nS, clustCtrX,clustCtrY,bestZ,whichClusters,whichBlocks)



h = figure;clf
set(h,'position',[1 1 1000 1000]);
ax = subplot(1,2,1,'parent',h);
hold(ax,'on')
colormap(ax,colormapRedBlue);
set(ax,'xdir','rev','zdir','rev','zlim',[1 51]);
%axis(ax,'vis3d')
threshold=.005;

 
for cc = whichClusters
    for bb = whichBlocks
        cf = load(nS.cellFileNameFcn(cc,bb),'rois');
        clustsThresholded = (sum(cf.rois,3)>threshold)*bestZ(cc,bb);
        clustsThresholded(clustsThresholded==0) = nan;
        
        frameSize = size(clustsThresholded);
        
        xStart = clustCtrX(cc,bb)-frameSize(2)/2+1/2;
        xEnd = clustCtrX(cc,bb)+frameSize(2)/2-1/2;
        yStart = clustCtrY(cc,bb)-frameSize(1)/2+1/2;
        yEnd = clustCtrY(cc,bb)+frameSize(1)/2-1/2;
        
        [xx, yy] = meshgrid(xStart:xEnd,yStart:yEnd);
        
        for addit=0, s = surf(ax,xx,yy,clustsThresholded+addit,round(clustsThresholded)); end
        alpha(s,.5);
        drawnow
    end
end
shading(ax,'flat')
grid(ax,'on')
xlabel(ax,'x'); ylabel(ax,'y'); zlabel(ax,'z')
set(ax,'color',[1 1 1]*.1)

