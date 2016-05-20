refPath = '/Users/tyler/Documents/Princeton/Rotations/tank/J115/2015-09-25__post_stack_002_AVERAGE.gif';
ref = imread(refPath,30);

%%
ref = cropStack(ref);
[refHeight refWidth] = size(ref);

blksz = 100;

blkStartY = 300;
blkStartX = 100;
blkEndY = blksz+blkStartY-1;
blkEndX = blksz+blkStartX-1;


img  = ref(101:400,101:400);
blk  = ref(blkStartY:blkEndY,blkStartX:blkEndX);

blkInRefLoc = false(size(ref));
blkInRefLoc(blkStartY:blkEndY,blkStartX:blkEndX) = true;
%%
angle = 10;
refRot = rotateAndSelectBlock(ref,true(size(ref)),angle);
noiseSD = 100;
noiseMu = 1;
movFrame = ref + cast(noiseSD .* randn(size(ref)) + noiseMu,'uint8');

[ blkRot] = rotateAndSelectBlock(movFrame,blkInRefLoc,angle);
corrMat   = computeBlockImageCorrs(blkRot, refRot,'double');


[yPeak xPeak] = ind2sub(size(corrMat),find(corrMat==max(corrMat(:))));
yPeakInRef = yPeak - blksz + 1;
xPeakInRef = xPeak - blksz + 1;

foundBlk   = refRot(yPeakInRef:yPeakInRef+blksz-1, xPeakInRef:xPeakInRef+blksz-1);

figure(1); clf
subplot(2,2,1)
imagesc(corrMat); axis image
title('correlation matrix');
subplot(2,2,2)
imshowpair(blkRot, foundBlk, 'falsecolor');
title('block differences');
axis image
subplot(2,2,3:4)
imshowpair(blkRot, foundBlk, 'montage');
title('blocks side by side')
set(gcf,'position',[100 100 1000 1000]); axis image