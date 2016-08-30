function  plotRoisInVol(nS, xyzrcoPeak,cc,bb,varargin)
p=inputParser;
addParameter(p,'roiThresh',.005);
addParameter(p,'whichRois',[]);
addParameter(p,'stackIm',[]);
parse(p,varargin{:});

stack = p.Results.stackIm;

zInRef = xyzrcoPeak(3,cc,bb);

% grab the found rois
cf  = load(nS.cellFileNameFcn(cc,bb),'rois');

cImSize = size(cf.rois);

% turn the center coords into indices of the whole block
xIndInRef = xyzrcoPeak(1,cc,bb) + (-cImSize(2)/2+1/2:cImSize(2)/2-1/2);
yIndInRef = xyzrcoPeak(2,cc,bb) +(-cImSize(1)/2+1/2:cImSize(1)/2-1/2);
[xx, yy]=meshgrid(xIndInRef,yIndInRef);

% figure out which pixels of the reference image to use and whether
% or not to pad the edges
refXInd = xIndInRef;
refXInd(xIndInRef<1)=[];
xLeftPad = sum(xIndInRef<1);
refYInd = yIndInRef;
refYInd(yIndInRef<1)=[];
yLeftPad = sum(yIndInRef<1);

% grab the relevant portion of the reference image
if ~isempty(stack)
    stackSliceBlockPadded = padarray(double((stack(round(refYInd),round(refXInd),zInRef))),[yLeftPad, xLeftPad], nan, 'pre');
end



nRois = size(cf.rois,3);

if ischar(p.Results.whichRois)  && strcmp(p.Results.whichRois,'all')
    whichRois = 1:nRois;
else
whichRois = p.Results.whichRois;
end



roiColorRange=linspace(0,1,nRois);
%roiColorRange = linspace(min(stackSliceBlockPadded(:)),max(stackSliceBlockPadded(:)), nRois);

smallZOffset = cc*.1;

if ~isempty(stack),
    hold on
    s=surf(xx,yy,zInRef*ones(size(xx)),normalizeToZeroOne(stackSliceBlockPadded))
    shading flat
    colormap bone
freezeColors

end
cmap = hsv(nRois);
colormap(cmap);
for rr = whichRois
    
    thisRoi = cf.rois(:,:,rr);
    thisRoi(thisRoi>p.Results.roiThresh)=1;
    thisRoi(thisRoi~=1)=nan;
    
    hold on
    
    roiS = surf(xx,yy,(zInRef)*ones(size(xx)),roiColorRange(rr).*thisRoi);%double(rois).*double(max(stackSliceBlockPadded(:))));
    alpha(roiS,.85);
    shading flat
    
    %roiXProfile = find(nansum(thisRoi,2)>0);
    %plot(clustsXZAx, roiXProfile, ones(size(roiXProfile)).*(zInRef+rand/10), 'color',cmap(rr,:));
end
%

%         %shading flat
%         %colormap(bone)


axis vis3d
axis ij
set(gca,'zdir','rev')








%
% h = figure;clf
% set(h,'position',[1 1 1000 1000]);
% ax = subplot(1,2,1,'parent',h);
% hold(ax,'on')
% colormap(ax,colormapRedBlue);
% set(ax,'xdir','rev','zdir','rev','zlim',[1 51]);
% %axis(ax,'vis3d')
% threshold=.005;
%
%
% for cc = whichClusters
%     for bb = whichBlocks
%         cf = load(nS.cellFileNameFcn(cc,bb),'rois');
%         clustsThresholded = (sum(cf.rois,3)>threshold)*bestZ(cc,bb);
%         clustsThresholded(clustsThresholded==0) = nan;
%
%         frameSize = size(clustsThresholded);
%
%         xStart = clustCtrX(cc,bb)-frameSize(2)/2+1/2;
%         xEnd = clustCtrX(cc,bb)+frameSize(2)/2-1/2;
%         yStart = clustCtrY(cc,bb)-frameSize(1)/2+1/2;
%         yEnd = clustCtrY(cc,bb)+frameSize(1)/2-1/2;
%
%         [xx, yy] = meshgrid(xStart:xEnd,yStart:yEnd);
%
%         for addit=0, s = surf(ax,xx,yy,clustsThresholded+addit,round(clustsThresholded)); end
%         alpha(s,.5);
%         drawnow
%     end
% end
% shading(ax,'flat')
% grid(ax,'on')
% xlabel(ax,'x'); ylabel(ax,'y'); zlabel(ax,'z')
% set(ax,'color',[1 1 1]*.1)

