refPath = '/Users/tyler/Documents/Princeton/Rotations/tank/J115/2015-09-25__post_stack_002_AVERAGE.gif';
ref = imread(refPath,30);

%%
ref = cropStack(ref);
[refHeight refWidth] = size(ref);

blksz = 100;

blkStartY = 375;
blkStartX = 200;
blkEndY = blksz+blkStartY-1;
blkEndX = blksz+blkStartX-1;


img  = ref(101:400,101:400);
blk  = ref(blkStartY:blkEndY,blkStartX:blkEndX);

blkInRefLoc = false(size(ref));
blkInRefLoc(blkStartY:blkEndY,blkStartX:blkEndX) = true;
%%
angle = 10;
refRot = imrotate(ref,angle,'bilinear');
noiseSD = 30;
noiseMu = 4;
movFrame = ref + cast(noiseSD .* randn(size(ref)) + noiseMu,'uint8');

[cm blkRot] = computeBlockImageCorrs(movFrame,blkInRefLoc,angle,refRot,'double');

[yPeak xPeak] = ind2sub(size(cm),find(cm==max(cm(:))));
yPeakInRef = yPeak - blksz + 1;
xPeakInRef = xPeak - blksz + 1;

foundBlk   = refRot(yPeakInRef:yPeakInRef+blksz-1, xPeakInRef:xPeakInRef+blksz-1);

figure(1); clf
imshowpair(blkRot, foundBlk, 'falsecolor');
set(gcf,'position',[100 100 1000 1000]); axis image