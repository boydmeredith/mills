
function plotRegRoiInVol(regRoi,clusterInd,cellDuplicateLists,...
    xyzrcoClusterPeaks,blockInd, stack)
% regRoi = regRoi1103;
% clusterInd=clusterInd1103;
% cellDuplicateLists=cellDuplicateLists1103;
% xyzrcoClusterPeaks=xyzrcoClusterPeaks1103;
% blockInd=blockInd1103;
% %%

ncolors = 20;
roiThresh=.001;

[stackH stackW stackD] = size(stack);
[xxStack yyStack] = meshgrid(1:stackW,1:stackH);



blockX = 100:600;
blockY = 100:600;
regRoiSz = size(regRoi);
roiSetInd = find(cellfun(@(x) length(x)>1, cellDuplicateLists));
cmap = hsv(ncolors);
[xx, yy]=meshgrid(1:length(blockX),1:length(blockY));
nRois = length(cellDuplicateLists);
roiColorRange=linspace(0,1,ncolors);

% subplot(1,2,1)
zToPlot = 29; surf(xxStack,yyStack,zToPlot.*ones(size(xxStack)),stack(:,:,zToPlot));
colormap bone
set(gca,'zdir','rev','ydir','rev','zlim',[25 45],...
            'ylim',[1 100],'xlim',[1 100], 'view',[0 90],'clim',[0 1]);

title('z = 29','fontsize',20)
freezeColors
hold on


axis image

axis([1 100 1 100])
shading flat

axis off

% % 
% % subplot(1,2,2)
% % zToPlot = 28; surf(xxStack,yyStack,zToPlot.*ones(size(xxStack)),stack(:,:,zToPlot));
% % colormap bone
% % set(gca,'zdir','rev','ydir','rev','zlim',[25 45],...
% %             'ylim',[1 100],'xlim',[1 100], 'view',[0 90],'clim',[0 1]);
% % freezeColors
% % hold on
% % 
% % title('z = 28','fontsize',20)
% % freezeColors
% % 
% % shading flat
% % 
% % axis image
% % 
% % axis off
% 
% 
% axis([1 100 1 100])
colormap(cmap);


pause()

for mm = 1:length(roiSetInd)
    
    
    thisRoiSetInd = roiSetInd(mm);
    roiInd = cellDuplicateLists{thisRoiSetInd};
    
    for rr = 1:length(roiInd)
        thisBlockNum = blockInd(roiInd(rr));
        thisClusterNum = clusterInd(roiInd(rr));

%         if thisClusterNum == 1
%             subplot(1,2,1)
%         elseif thisClusterNum==2
%             subplot(1,2,2)
%         end
%         set(gca,'zdir','rev','ydir','rev','zlim',[25 45],...
%             'ylim',[1 100],'xlim',[1 100], 'view',[0 90]);
        
        thisClusterZ = xyzrcoClusterPeaks(3,thisClusterNum,thisBlockNum);
        
        thisRoi = regRoi(blockY,blockX,roiInd(rr));
        thisRoi(thisRoi<roiThresh) = nan;
        thisRoi(thisRoi>=roiThresh) = 1;
        surf(xx,yy,(thisClusterZ-.1*thisClusterNum).*ones(size(xx)),(roiColorRange(mod(mm,ncolors)+1)).*thisRoi);
        shading flat
        pause()
        %alpha(.7)
        
    end
    
    hold on
end