function [ballStickFig] = xyzrcBallStickFig(xyzrcPeak, thisFrameNo, ballStickFig, stackDim)


nBlocks = size(xyzrcPeak,2);
nFrames = size(xyzrcPeak,3);
frameString = sprintf('frame: %03i/%03i',thisFrameNo,nFrames);

if isempty(ballStickFig), ballStickFig = figure; end;
if isempty(stackDim), 
    xlimits = [-100, 600]; 
    ylimits = [-100, 600]; 
    zlimits = [min(reshape(xyzrcPeak(3,:),[],1)), max(reshape(xyzrcPeak(3,:),[],1))];
else
    xlimits = [-100 stackDim.width*(1+1/sqrt(nBlocks))];
    ylimits = [-100 stackDim.height*(1+1/sqrt(nBlocks))]; 
    zlimits = [1 stackDim.depth];
end


% ball and stick image plot circles for the block centers with size
%   proportional to rotation index and colored by their correlation value
thisFramePeaks = xyzrcPeak(:,:,thisFrameNo);
thisFramePeaksGridX = reshape(thisFramePeaks(1,:), sqrt(nBlocks), sqrt(nBlocks));
thisFramePeaksGridY = reshape(thisFramePeaks(2,:), sqrt(nBlocks), sqrt(nBlocks));
thisFramePeaksGridZ = reshape(thisFramePeaks(3,:), sqrt(nBlocks), sqrt(nBlocks));

figure(clf(ballStickFig))
ballStickAx = axes();

markersizes = eps+normalizeToZeroOne(reshape(thisFramePeaks(5,:),[],1)) * 8e2;%5e3*(reshape(thisFramePeaks(4,:),[],1)+abs(min(reshape(thisFramePeaks(4,:),[],1)))) + eps;

scatter3(ballStickAx, reshape(thisFramePeaks(1,:),[],1),reshape(thisFramePeaks(2,:),[],1),...
    reshape(thisFramePeaks(3,:),[],1),...
    markersizes,...%max(reshape(thisFramePeaks(4,:),[],1),1),...
    reshape(thisFramePeaks(4,:),[],1) , 'filled','linewidth',2,'markeredgecolor',[.3 .3 .3]);


hold(ballStickAx, 'on');
% connect the circles with lines
for bb = 1:sqrt(nBlocks)
    plot3(ballStickAx, thisFramePeaksGridX(bb,:),thisFramePeaksGridY(bb,:),...
        thisFramePeaksGridZ(bb,:),'-','color',[.4 .4 .5]);
    plot3(ballStickAx, thisFramePeaksGridX(:,bb),thisFramePeaksGridY(:,bb),...
        thisFramePeaksGridZ(:,bb),'-','color',[.4 .4 .5]);
end
set(ballStickAx,'xdir','rev','xlim',xlimits,...
    'ylim',ylimits,...
    'zdir','reverse','zlim',zlimits,...
    'xaxislocation','origin','yaxislocation','origin',...
    'CameraPosition',[37.5032 2.7521e+03 -43.6450]);
axis(ballStickAx,'square');
set(ballStickFig,'Position',[560         227        1029         721]);
xlabel(ballStickAx,'x'); ylabel(ballStickAx,'y'); zlabel(ballStickAx,'z');
colormap(ballStickAx, colormapRedBlue); colorbar(ballStickAx);

title(ballStickAx, frameString);

set(gca,'color',[.8 .8 .8]);


end