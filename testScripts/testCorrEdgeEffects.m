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
figure(2); clf
h = makeSubplots(2,3,3,.1,.2,[0 0 1 1]);
axes(h(1));
imagesc(ref); axis image; colormap(bone); freezeColors;
title('reference');
axes(h(2));
imagesc(img); axis image; colormap(bone); freezeColors;
title('image');

axes(h(3));
title('block');
imagesc(blk); axis image;
colormap(bone); freezeColors;

axes(h(5));
imagesc(normxcorr2(blk,img)); axis image
colormap(redbluecmap); freezeColors;
caxis([-1 1]);
title('normxcorr2(blk,img)');

axes(h(6));
imagesc(xcorr2(img,blk)); axis image
colormap(redbluecmap); freezeColors;
title('xcorr2(img,blk)');

% iterate through block positions and take correlation of overlap
refPadded = padarray(ref, [blksz-1, blksz-1], -Inf);
refLocPadded = padarray(true(size(ref)), [blksz-1, blksz-1], false);

refToImg = -Inf*ones(size(ref));
refToImg(101:400,101:400) = ref(101:400,101:400);

imgPadded = padarray(refToImg, [blksz-1, blksz-1], -Inf);;


imgLocPadded = true(size(imgPadded));
imgLocPadded(imgPadded == -Inf) = false;

%%

axes(h(8));
imagesc(refLocPadded); axis image; colormap(bone); freezeColors;
corrMat = zeros(refHeight+blksz-1,refWidth+blksz-1);

for sx = 1:(refWidth+blksz-1) %(blkStartX+blksz-1):(blkStartX+blksz-1);
    for sy = 1:(refHeight+blksz-1) %(blkStartY+blksz-1):(blkStartY+blksz-1)%;
        fprintf('x=%i, y=%i...  \t',sx,sy);
        
        blkLocPadded = false(size(refPadded));
        blkLocPadded(sy:(sy+blksz-1),sx:(sx+blksz-1)) = true;
        
        refOverlap = refLocPadded & blkLocPadded;
        [rY rX]    = getBlockIx(refOverlap);
        
        blkRefOverlap = refOverlap(sy:(sy+blksz-1),sx:(sx+blksz-1));
        [bRY bRX]    = getBlockIx(blkRefOverlap);
        
        imgOverlap = imgLocPadded & blkLocPadded;
        [iY iX]    = getBlockIx(imgOverlap);
        
        blkImgOverlap = imgOverlap(sy:(sy+blksz-1),sx:(sx+blksz-1));
        [bIY bIX]    = getBlockIx(blkImgOverlap);
        
        
        refChunk = refPadded(rY,rX);
        imgChunk = imgPadded(iY, iX);
        blkRefChunk = blk(bRY,bRX);
        blkImgChunk = blk(bIY,bIX);
        %         if sx == sy & sum(blkOverlap(:)) > 1
        %             axes(h(8)); cla
        %             imagesc(refOverlap); axis image; colormap(bone); freezeColors;
        %             drawnow
        %             axes(h(6)); cla
        %             imagesc(corrMat);
        %             title('corr on overlap');
        %             axis image
        %             colormap(redbluecmap); freezeColors;
        % %             imagesc(blkOverlap); axis image; colormap(bone); freezeColors;
        %             drawnow
        %         end
        
        
        % corrMat(sx,sy) = corr(cast(blkChunk(:),'double'),cast(refChunk(:),'double'))
        if sx == (blkStartX+blksz-1) & sy == (blkStartY+blksz-1)
            axes(h(4));
            imagesc(refLocPadded+blkLocPadded+imgLocPadded)
            axis image
        end
        if ~isempty(blkImgChunk)
            corrMat(sy,sx) = corr(cast(blkImgChunk(:),'double'),cast(imgChunk(:),'double'));
        end
    end
end

%%

[yPeak xPeak] = ind2sub(size(corrMat),find(corrMat==max(corrMat(:))));

axes(h(6)); cla
imagesc(flipud(corrMat(blksz-1:end-blksz+1,blksz-1:end-blksz+1)));
hold on
title('corr on overlap');
scatter(xPeak,yPeak,'cx');
axis image

colormap(redbluecmap);
xlim([1 400])
ylim([1 400])

axes(h(4)); colorbar
axes(h(6)); colorbar

