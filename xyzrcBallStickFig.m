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
hold(ballStickAx, 'on');

if size(xyzrcPeak,1) == 5
    corrValsZeroOne = zeros(squeeze(size(xyzrcPeak(5,:,:))));
    corrValsZeroOne(xyzrcPeak(5,:)~=0) = squeeze(normalizeToZeroOne(xyzrcPeak(5,xyzrcPeak(5,:)~=0)));
    allmarkersizes = squeeze(100+corrValsZeroOne .* 7e2);%5e3*(reshape(thisFramePeaks(4,:),[],1)+abs(min(reshape(thisFramePeaks(4,:),[],1)))) + eps;
    markersizes = allmarkersizes(:,thisFrameNo);
    bubsizes = [min(allmarkersizes(xyzrcPeak(5,:)~=0)) quantile(allmarkersizes(xyzrcPeak(5,:)~=0),[0.5]) max(allmarkersizes(xyzrcPeak(5,:)~=0))];
    bubsizelabels = [min(xyzrcPeak(5,xyzrcPeak(5,:)~=0)) quantile(xyzrcPeak(5,xyzrcPeak(5,:)~=0),[0.5]) max(xyzrcPeak(5,xyzrcPeak(5,:)~=0))];
    legentry=cell(size(bubsizes));
    for ind = 1:numel(bubsizes)
       bubleg(ind) = plot(ballStickAx, 0,0,'ro','markersize',sqrt(bubsizes(ind)),'MarkerFaceColor','red');
       set(bubleg(ind),'visible','off')
       legentry{ind} = num2str(round(bubsizelabels(ind),2));
    end
    hLegend = legend(ballStickAx,legentry);
    set(hLegend,'box','off');
    axl = axes('Parent',hLegend.Parent, 'Units',hLegend.Units, 'Position',hLegend.Position, ...
                       'XTick',[] ,'YTick',[], 'Color','none', 'YColor','none', 'XColor','none', 'HandleVisibility','off', 'HitTest','off');
    title(axl,'correlation');
    climits = [min(xyzrcPeak(4,xyzrcPeak(5,:)~=0)) max(xyzrcPeak(4,xyzrcPeak(5,:)~=0))];
else
    markersizes = 500;
    climits = [min(xyzrcPeak(4,:)) max(xyzrcPeak(4,:))];
end


scatter3(ballStickAx, thisFramePeaks(1,:),thisFramePeaks(2,:),...
    thisFramePeaks(3,:),...
    markersizes,...%max(reshape(thisFramePeaks(4,:),[],1),1),...
    thisFramePeaks(4,:) , 'filled','linewidth',2,'markeredgecolor',[.3 .3 .3]);


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
    'CameraPosition',[37.5032 2.7521e+03 -43.6450]);
axis(ballStickAx,'square');
set(ballStickFig,'Position',[560         227        1029         721]);
xlabel(ballStickAx,'x'); ylabel(ballStickAx,'y'); zlabel(ballStickAx,'z');
colormap(ballStickAx, hot); c= colorbar(ballStickAx);
c.Label.String = 'Rotation Angle';






title(ballStickAx, frameString);

set(gca,'color',[.8 .8 .8],'clim',climits);


end