function [ballStickFig] = xyzrcBallStickFig(xyzrcoPeak, thisFrameNo, ballStickFig, stackDim, mByN)


nBlocks = size(xyzrcoPeak,2);
nFrames = size(xyzrcoPeak,3);
frameString = sprintf('frame: %03i/%03i',thisFrameNo,nFrames);

if isempty(ballStickFig), ballStickFig = figure; end;
if isempty(stackDim), 
    xlimits = [-100, 600]; 
    ylimits = [-100, 600]; 
    zlimits = [min(reshape(xyzrcoPeak(3,:),[],1)), max(reshape(xyzrcoPeak(3,:),[],1))];
else
    xlimits = [-100 stackDim.width*(1+1/mByN(2))];
    ylimits = [-100 stackDim.height*(1+1/mByN(1))]; 
    zlimits = [1 stackDim.depth];
end


% ball and stick image plot circles for the block centers with size
%   proportional to rotation index and colored by their correlation value
thisFramePeaks = xyzrcoPeak(:,:,thisFrameNo);
thisFramePeaksGridX = reshape(thisFramePeaks(1,:), mByN(1), mByN(2));
thisFramePeaksGridY = reshape(thisFramePeaks(2,:), mByN(1), mByN(2));
thisFramePeaksGridZ = reshape(thisFramePeaks(3,:), mByN(1), mByN(2));

clf(ballStickFig);
ballStickAx = axes('parent',ballStickFig);
hold(ballStickAx, 'on');

if size(xyzrcoPeak,1) >= 5
    corrValsZeroOne = zeros(squeeze(size(xyzrcoPeak(5,:,:))));
    corrValsZeroOne(xyzrcoPeak(5,:)~=0) = squeeze(normalizeToZeroOne(xyzrcoPeak(5,xyzrcoPeak(5,:)~=0)));
    allmarkersizes = squeeze(100+corrValsZeroOne .* 7e2);%5e3*(reshape(thisFramePeaks(4,:),[],1)+abs(min(reshape(thisFramePeaks(4,:),[],1)))) + eps;
    markersizes = allmarkersizes(:,thisFrameNo);
    bubsizes = [min(allmarkersizes(xyzrcoPeak(5,:)~=0)) quantile(allmarkersizes(xyzrcoPeak(5,:)~=0),[0.5]) max(allmarkersizes(xyzrcoPeak(5,:)~=0))];
    bubsizelabels = [min(xyzrcoPeak(5,xyzrcoPeak(5,:)~=0)) quantile(xyzrcoPeak(5,xyzrcoPeak(5,:)~=0),[0.5]) max(xyzrcoPeak(5,xyzrcoPeak(5,:)~=0))];
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
    climits = [min(xyzrcoPeak(4,xyzrcoPeak(5,:)~=0)) max(xyzrcoPeak(4,xyzrcoPeak(5,:)~=0))];
    

else
    markersizes = 500;
    climits = [min(xyzrcoPeak(4,:)) max(xyzrcoPeak(4,:))];
end


%surf(thisFramePeaksGridX,thisFramePeaksGridY,thisFramePeaksGridZ,'facecolor','blue','parent',ballStickAx);
%alpha(.3)

scatter3(ballStickAx, thisFramePeaks(1,:),thisFramePeaks(2,:),...
    thisFramePeaks(3,:),...
    markersizes,...%max(reshape(thisFramePeaks(4,:),[],1),1),...
    thisFramePeaks(4,:) , 'filled','linewidth',2,'markeredgecolor',[.8 .8 .8]);


% add dark edges to circles referring to outliers
if size(xyzrcoPeak,1) == 6,
    outliers = find(thisFramePeaks(6,:));
    scatter3(ballStickAx, thisFramePeaks(1,outliers),thisFramePeaks(2,outliers),...
    thisFramePeaks(3,outliers),...
    markersizes(outliers),...%max(reshape(thisFramePeaks(4,:),[],1),1),...
    thisFramePeaks(4,outliers) , 'filled','linewidth',2,'markeredgecolor','k');
end


% connect the circles with lines
plot3(ballStickAx, thisFramePeaksGridX,thisFramePeaksGridY,...
        thisFramePeaksGridZ,'-','color',[.4 .4 .6]);
plot3(ballStickAx, thisFramePeaksGridX',thisFramePeaksGridY',...
        thisFramePeaksGridZ','-','color',[.4 .4 .6]);

set(ballStickAx,'xdir','rev','xlim',xlimits,...
    'ylim',ylimits,...
    'zdir','reverse','zlim',zlimits,...
    'CameraPosition',[37.5032 2.7521e+03 -43.6450],...
    'clim',climits);
axis(ballStickAx,'square');
set(ballStickFig,'Position',[560         227        1029         721]);
xlabel(ballStickAx,'x'); ylabel(ballStickAx,'y'); zlabel(ballStickAx,'z');
colormap(ballStickAx, hot); c= colorbar(ballStickAx);
c.Label.String = 'Rotation Angle';



grid(ballStickAx,'on')


title(ballStickAx, frameString);



end