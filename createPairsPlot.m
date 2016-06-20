

function [pairsPlot] = createPairsPlot(figNum, blockSt, title, cLimits)
% X-Y
[~, peakInd] = max(blockSt.corrValsToSave(:));
clf(figNum);

if ~strcmp(class(blockSt.corrValsToSave),'double') 
    blockSt.corrValsToSave = double(blockSt.corrValsToSave)/ double(intmax(class(blockSt.corrValsToSave)));
end

refHeight = 512;
refWidth = 512;

pairsPlot = makeSubplots(figNum,3,2,.1,.1,[0 0 1 1]);
axes(pairsPlot(1));

blockSt.indZRPeak = find(blockSt.indZ==blockSt.indZ(peakInd) & blockSt.indR ==blockSt.indR(peakInd));
imagesc(sparse(double(blockSt.indY(blockSt.indZRPeak)),double(blockSt.indX(blockSt.indZRPeak)),  ...
    blockSt.corrValsToSave(blockSt.indZRPeak), refHeight,refWidth) );
ylabel('Y');
xlabel('X');
axis image; colorbar; colormap(colormapRedBlue); caxis(cLimits);

% X-R
axes(pairsPlot(2));
blockSt.indZYPeak = find(blockSt.indZ==blockSt.indZ(peakInd) & blockSt.indY ==blockSt.indY(peakInd));
imagesc(sparse(double(blockSt.indR(blockSt.indZYPeak)), double(blockSt.indX(blockSt.indZYPeak)), ...
    blockSt.corrValsToSave(blockSt.indZYPeak)));
ylabel('R');
xlabel('X');
axis image; colorbar; colormap(colormapRedBlue); caxis(cLimits);

% plot title
text(mean(get(gca,'xlim')), -15,title,'fontsize',25);

% X-Z
axes(pairsPlot(3));
blockSt.indYRPeak = find(blockSt.indY==blockSt.indY(peakInd) & blockSt.indR ==blockSt.indR(peakInd));
imagesc(sparse(double(blockSt.indZ(blockSt.indYRPeak)), double(blockSt.indX(blockSt.indYRPeak)), ...
    blockSt.corrValsToSave(blockSt.indYRPeak)));
ylabel('R');
xlabel('X');
axis image; colorbar; colormap(colormapRedBlue); caxis(cLimits);

% Y-R
axes(pairsPlot(4));
blockSt.indZXPeak = find(blockSt.indZ==blockSt.indZ(peakInd) & blockSt.indX ==blockSt.indX(peakInd));
imagesc(sparse(double(blockSt.indR(blockSt.indZXPeak)), double(blockSt.indY(blockSt.indZXPeak)),  ...
    blockSt.corrValsToSave(blockSt.indZXPeak)));
ylabel('R');
xlabel('Y');
axis image; colorbar; colormap(colormapRedBlue); caxis(cLimits);

% Y-Z
axes(pairsPlot(5));
blockSt.indXRPeak = find(blockSt.indX==blockSt.indX(peakInd) & blockSt.indR ==blockSt.indR(peakInd));
imagesc(sparse(double(blockSt.indZ(blockSt.indXRPeak)), double(blockSt.indY(blockSt.indXRPeak)), ...
    blockSt.corrValsToSave(blockSt.indXRPeak)));
xlabel('Y');
ylabel('Z');
axis image; colorbar; colormap(colormapRedBlue); caxis(cLimits);

% Z-R
axes(pairsPlot(6));
blockSt.indXYPeak = find(blockSt.indX==blockSt.indX(peakInd) & blockSt.indY ==blockSt.indY(peakInd));
imagesc(sparse(double(blockSt.indR(blockSt.indXYPeak)), double(blockSt.indZ(blockSt.indXYPeak)), ...
    blockSt.corrValsToSave(blockSt.indXYPeak)));
ylabel('R');
xlabel('Z');
axis image; colorbar; colormap(colormapRedBlue); caxis(cLimits);
end