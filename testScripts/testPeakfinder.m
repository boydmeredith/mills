nxc2M = [];
nxc2MNbrhd = [];
bestBlockRot = rotateAndSelectBlock(movieFrame, thisBlockLoc, params.rotAngleFromInd(indR(peakInd)));

%%
NXYTOKEEP = 10000;


for zz = 1:length(zRangeToKeep)
    stackSlice = normalizeToZeroOne(stack(:,:,zRangeToKeep(zz)));
    nxc2MNbrhd(:,:,:,zz) = computeBlockImageCorrs(bestBlockRot, ...
        stackSlice, nbrhdInf, params.minCorrOverlap, 'double');
    nxc2M(:,:,:,zz) = normxcorr2(bestBlockRot, stackSlice);
    
    
    [yIx, xIx, thisCorr] = find(nxc2MNbrhd(:,:,:,zz));
    [~, thisCorrSortIx] = sort(thisCorr,'descend');
    thisCorrXYToKeepIx  = thisCorrSortIx(1:NXYTOKEEP);
    
    figure(7); clf;
    subplot(2,3,4)
    imagesc(bestBlockRot); axis image
    subplot(2,3,5:6)
    imagesc(stackSlice); axis image
    hold on
    text(xyzrPrior(1,thisBlockNo),...
        xyzrPrior(2,thisBlockNo),'P','fontsize',20)
    
    subplot(2,3,1)
    imagesc(nxc2M(:,:,:,zz)); title(num2str(zRangeToKeep(zz)));
    caxis([min(nxc2M(:)) max(nxc2M(:))]);
    colormap(colormapRedBlue);
    [yPeak, xPeak] = find(max(reshape(nxc2M(:,:,:,zz),[],1))==nxc2M(:,:,:,zz))
    hold on
    scatter(xPeak, yPeak, 200,'ko')
    text(xyzrPrior(1,thisBlockNo)+bInf.width/2-1/2,...
        xyzrPrior(2,thisBlockNo)+bInf.height/2-1/2,'P','fontsize',20)
    
    subplot(2,3,2)
    imagesc(nxc2MNbrhd(:,:,:,zz)); title(num2str(zRangeToKeep(zz)));
    caxis([min(nxc2M(:)) max(nxc2M(:))]);
    colormap(colormapRedBlue);
    [yPeakN, xPeakN] = find(max(reshape(nxc2MNbrhd(:,:,:,zz),[],1))==nxc2MNbrhd(:,:,:,zz))
    hold on
    scatter(xPeakN, yPeakN, 200,'ko')
    text(xyzrPrior(1,thisBlockNo)+bInf.width/2-1/2,...
        xyzrPrior(2,thisBlockNo)+bInf.height/2-1/2,'P','fontsize',20)
    
    subplot(2,3,3)
    imagesc(sparse(yIx(thisCorrXYToKeepIx), xIx(thisCorrXYToKeepIx),...
        thisCorr(thisCorrXYToKeepIx), size(nxc2M(:,:,:,zz),1),size(nxc2M(:,:,:,zz),2)));
    caxis([min(nxc2M(:)) max(nxc2M(:))]);
    colormap(colormapRedBlue);
    hold on
    scatter(xIx(thisCorrXYToKeepIx(1)), yIx(thisCorrXYToKeepIx(1)), 200,'ko')
    text(xyzrPrior(1,thisBlockNo)+bInf.width/2-1/2,...
        xyzrPrior(2,thisBlockNo)+bInf.height/2-1/2,'P','fontsize',20)
    pause()
    
    
end

is = imageSeries(nxc2MNbrhd);
is.plotImagesSequentially;
