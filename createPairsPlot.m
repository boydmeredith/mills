function [pairsPlot] = createPairsPlot(subj, movieDate, blockNo, frameNo)
% function [pairsPlot] = createPairsPlot(subj, movieDate, blockNo, frameNo)
%
% Plots correlation values around the peak for 2 and 1 dimension
% combinations of x,y,z,r
% 
% input:
% subj      - e.g. 'J114'
% movieDate - e.g. '2015-09-17'
% blockNo   - which block to plot correlations for
% frameNo   - which frame of the movie to plot correlations for

figTitle = sprintf('%s %s block: %03d frame %03d', subj, movieDate, blockNo, frameNo);

blockPath = fullfile(jlgDataDir, subj, movieDate,...
    sprintf('referenceLocalization/frame%03d/block%03d.mat', frameNo, blockNo));

blockSt = load(blockPath);

refinfo = imfinfo(fullfile(jlgDataDi, blockSt.stackPath));

refDepth = length(refinfo);
refHeight = refinfo(1).Height;
refWidth = refinfo(1).Width;



[~, peakInd] = max(blockSt.corrValsToSave(:));
h = figure; 
clf(h);
set(h, 'Position', [90  26 1401 922]);

if ~strcmp(class(blockSt.corrValsToSave),'double') 
    blockSt.corrValsToSave = double(blockSt.corrValsToSave)/ double(intmax(class(blockSt.corrValsToSave)));
end

refHeight = 512;
refWidth = 512;

pairsPlot = makeSubplots(h,3,2,.2,.35,[0 .3 1 .7]);
linesPlot = makeSubplots(h,4,1,.2,.2,[.05 0 .9 .3]);

xPeak = blockSt.indX(peakInd);
yPeak = blockSt.indY(peakInd);
zPeak = blockSt.indZ(peakInd);
rPeak = blockSt.indR(peakInd);

indXPeak = find(blockSt.indX == xPeak);
indYPeak = find(blockSt.indY == yPeak);
indZPeak = find(blockSt.indZ == zPeak);
indRPeak = find(blockSt.indR == rPeak);

indZRPeak = intersect(indZPeak,indRPeak);
indZYPeak = intersect(indZPeak,indYPeak);
indYRPeak = intersect(indYPeak,indRPeak);
indZXPeak = intersect(indZPeak,indXPeak);
indXYPeak = intersect(indYPeak,indXPeak);
indXRPeak = intersect(indXPeak,indRPeak);


indZRYPeak = intersect(indZRPeak,indYPeak);
indZRXPeak = intersect(indZRPeak,indXPeak);
indYRXPeak = intersect(indYRPeak,indXPeak);
indZYXPeak = intersect(indZYPeak,indXPeak);


[~, s] = sort(blockSt.indX(indZRYPeak));
plot(linesPlot(1,1), blockSt.indX(indZRYPeak(s)), blockSt.corrValsToSave(indZRYPeak(s)),'o-',...
    repmat(xPeak,1,2),quantile(blockSt.corrValsToSave(indZRYPeak(s)),[0 1]),'r'); 
ylabel(linesPlot(1,1),'correlation')
xlabel(linesPlot(1,1),'X')

[~, s] = sort(blockSt.indY(indZRXPeak));
plot(linesPlot(2,1), blockSt.indY(indZRXPeak(s)), blockSt.corrValsToSave(indZRXPeak(s)),'o-',...
    repmat(yPeak,1,2),quantile(blockSt.corrValsToSave(indZRXPeak(s)),[0 1]),'r');
xlabel(linesPlot(2,1),'Y');

[~, s] = sort(blockSt.indZ(indYRXPeak));
plot(linesPlot(3,1), blockSt.indZ(indYRXPeak(s)), blockSt.corrValsToSave(indYRXPeak(s)),'o-',...
    repmat(zPeak,1,2),quantile(blockSt.corrValsToSave(indYRXPeak(s)),[0 1]),'r');
xlabel(linesPlot(3,1),'Z');

[~, s] = sort(blockSt.indR(indZYXPeak));
plot(linesPlot(4,1), blockSt.rotAngleFromInd(blockSt.indR(indZYXPeak(s))), ...
    blockSt.corrValsToSave(indZYXPeak(s)),'o-',...
    repmat(blockSt.rotAngleFromInd(rPeak),1,2),quantile(blockSt.corrValsToSave(indZYXPeak(s)),[0 1]),'r'); 
xlabel(linesPlot(4,1),'R');



% X-Y
xyAx = pairsPlot(1);
imagesc(sparse(double(blockSt.indY(indZRPeak)),double(blockSt.indX(indZRPeak)),  ...
    blockSt.corrValsToSave(indZRPeak), refHeight,refWidth),'parent',xyAx );
ylabel(xyAx,'Y');
xlabel(xyAx,'X');
axis(xyAx, 'square'); colorbar(xyAx); 
axis(xyAx,double([min(blockSt.indX(indZRPeak)) max(blockSt.indX(indZRPeak)) ...
    min(blockSt.indY(indZRPeak)) max(blockSt.indY(indZRPeak))])+[-10 10 -10 10]);
hold(xyAx,'on'); scatter(xyAx,xPeak,yPeak,'x');

% X-Z
xzAx = pairsPlot(2);
imagesc(sparse(double(blockSt.indZ(indYRPeak)), double(blockSt.indX(indYRPeak)), ...
    blockSt.corrValsToSave(indYRPeak), refDepth, refWidth),'parent',xzAx);
ylabel(xzAx,'Z');
xlabel(xzAx,'X');
axis(xzAx, 'square'); colorbar(xzAx); 
axis(xzAx,[-10 10 -10 10]+double([min(blockSt.indX(indYRPeak)) max(blockSt.indX(indYRPeak)) ...
    min(blockSt.indZ(indYRPeak)) max(blockSt.indZ(indYRPeak))])) 
hold(xzAx,'on'); scatter(xzAx,xPeak,zPeak,'x');

% plot title
title(pairsPlot(2), figTitle);

% X-R
xrAx = pairsPlot(4);
imagesc(sparse(double(blockSt.indR(indZYPeak)), double(blockSt.indX(indZYPeak)), ...
    blockSt.corrValsToSave(indZYPeak), length(blockSt.rotAngleFromInd), refWidth),'parent',xrAx);
ylabel(xrAx,'R');
xlabel(xrAx,'X');
axis(xrAx, 'square'); colorbar(xrAx); 
axis(xrAx,[-10 10 -10 10]+double([min(blockSt.indX(indZYPeak)) max(blockSt.indX(indZYPeak))...
    min(blockSt.indR(indZYPeak)) max(blockSt.indR(indZYPeak))])) 
hold(xrAx,'on'); scatter(xrAx,xPeak,rPeak,'x');





% Y-R
yrAx = pairsPlot(5);
imagesc(sparse(double(blockSt.indR(indZXPeak)), double(blockSt.indY(indZXPeak)),  ...
    blockSt.corrValsToSave(indZXPeak),length(blockSt.rotAngleFromInd), refHeight),'parent',yrAx);
ylabel(yrAx,'R');
xlabel(yrAx,'Y');
axis(yrAx, 'square'); colorbar(yrAx); 
axis(yrAx,[-10 10 -10 10]+double([min(blockSt.indY(indZXPeak)) max(blockSt.indY(indZXPeak)) ...
    min(blockSt.indR(indZXPeak)) max(blockSt.indR(indZXPeak))])) 
hold(yrAx,'on'); scatter(yrAx,yPeak,rPeak,'x');

% Y-Z
yzAx = pairsPlot(3)
imagesc(sparse(double(blockSt.indZ(indXRPeak)), double(blockSt.indY(indXRPeak)), ...
    blockSt.corrValsToSave(indXRPeak), refDepth, refHeight),'parent',yzAx);
xlabel(yzAx,'Y');
ylabel(yzAx,'Z');
axis(yzAx, 'square'); colorbar(yzAx); 
axis(yzAx,[-10 10 -10 10]+double([min(blockSt.indY(indXRPeak)) max(blockSt.indY(indXRPeak)) ...
    min(blockSt.indZ(indXRPeak)) max(blockSt.indZ(indXRPeak))])) 
hold(yzAx,'on'); scatter(yzAx,yPeak,zPeak,'x');

% Z-R
zrAx = pairsPlot(6);
imagesc(sparse(double(blockSt.indR(indXYPeak)), double(blockSt.indZ(indXYPeak)), ...
    blockSt.corrValsToSave(indXYPeak), length(blockSt.rotAngleFromInd), refDepth),'parent',zrAx);
hold(zrAx,'on')
axis(zrAx,[-10 10 -10 10]+double([min(blockSt.indZ(indXYPeak)) max(blockSt.indZ(indXYPeak)) ...
    min(blockSt.indR(indXYPeak)) max(blockSt.indR(indXYPeak))])) 
hold(zrAx,'on'); scatter(zrAx,zPeak,rPeak,'x');



ylabel(zrAx,'R');
xlabel(zrAx,'Z');
axis(zrAx, 'square'); colorbar(zrAx); 
end