climits = [min(xyzrcPeak(4,:)) max(xyzrcPeak(4,:))];
im = [], map = [];
h = figure(1);
clf(h)
set(h,'paperpositionmode','auto');
for ff = 1:63
h = xyzrcBallStickFig(xyzrcPeak,ff,h,[]);
set(h,'Position',[10 10 1000 1000])
set(findobj(h,'type','axes'),'CameraPosition', [250 250 -367.2114],...
    'View',[-179.6000 90], 'clim', climits)
[im, map] = createGif(h, ff, 63, im, map, 'toothDecay.gif');
end