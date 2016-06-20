climits = [min(xyzrcPeak(4,:)) max(xyzrcPeak(4,:))];
im = [], map = [];
birdsEye = true;
if birdsEye
    viewPoint = [-179.6000 90];
    fname = 'ballStickBirdsEye.gif';
else
    viewPoint = [166 12];
    fname = 'ballStickSide.gif';
end
stackDim.width = 512;
stackDim.height = 512;
stackDim.depth = 51;
for ff = 1:70
    h= figure;
    set(h,'Position',[10 10 1000 1000],'paperpositionmode','auto')

    h = xyzrcBallStickFig(xyzrcPeak,ff,h,stackDim);
    set(findobj(h,'type','axes'),'CameraPosition', [200 200 -300.2114],...
        'View',viewPoint, 'clim', climits)
    [im, map] = createGif(h, ff, 63, im, map, fname);
    close(h)
    pause(1)
    
end