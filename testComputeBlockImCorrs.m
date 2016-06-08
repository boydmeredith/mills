peppers = rgb2gray(imread('peppers.png'));
yStart = 1;
xStart = 1;
bWidth  = 70;
bHeight = 70;
block   = peppers(yStart:yStart+bHeight-1, xStart:xStart+bWidth-1);

corrMat = normxcorr2(block,peppers);
figure(1); clf
subplot(1,2,1)
imagesc(corrMat); axis square
colormap(redbluecmap)
colorbar
subplot(1,2,2)
colorbar

[yPeak, xPeak] = find(corrMat == max(corrMat(:)));
yUL = yPeak - size(block,1) + 1;
xUL = xPeak - size(block,2) + 1;
match  = peppers(yUL:yUL+size(block,1)-1, xUL:xUL+size(block,2)-1);
mean(block(:)-match(:))
subplot(1,3,1); imagesc(block);
subplot(1,3,2); imagesc(match);
subplot(1,3,3); imagesc(block - match);


%
figure(2); clf
nbrhdInf.yCtr = yUL + (size(block,1)-1)/2 - 10;
nbrhdInf.xCtr = xUL + (size(block,2)-1)/2 - 5;

yUL = nbrhdInf.yCtr - (size(block,1)-1)/2;
xUL = nbrhdInf.xCtr - (size(block,2)-1)/2;

% match  = peppers(yUL:yUL+size(onion,1)-1, xUL:xUL+size(onion,2)-1);
% mean(onion(:)-match(:))
%   
nbrhdInf.yMargin = 80;
nbrhdInf.xMargin = 70;
corrMatNbrhd = computeBlockImageCorrs(block,peppers,nbrhdInf,.5,'double');
subplot(1,3,1)
imagesc(corrMat); title('full corr mat'); axis square
subplot(1,3,2)
imagesc(corrMatNbrhd); title('nbrhd corr mat'); axis square
subplot(1,3,3); 
imagesc(corrMat-corrMatNbrhd); title('difference'); caxis([-.001 .001]); axis square
colormap(redbluecmap)

[realYPeak, realXPeak] = find(corrMat == max(corrMat(:)));
[myYPeak, myXPeak] = find(corrMatNbrhd == max(corrMatNbrhd(:)));
isequal([realYPeak, realXPeak], [myYPeak, myXPeak])




%%
blockWidth = size(block,2);
blockHeight = size(block,1);
overlap = .9;
minOverlap = blockHeight * blockWidth * overlap ; 
if overlap > 1
    minOverlap = overlap;
end
stackWidth = size(peppers,2);
stackHeight = size(peppers,1);

corrMatWidth = size(corrMat,2);
corrMatHeight = size(corrMat,1);

C = min(repmat(min(1:corrMatWidth,corrMatWidth:-1:1),corrMatHeight,1),   blockWidth);
R = min(repmat(min(1:corrMatHeight, corrMatHeight:-1:1)',1,corrMatWidth),blockHeight);

figure(3); clf
subplot(1,2,1); imagesc(corrMat);
corrMatOverlap = corrMat;
corrMatOverlap(R.*C < minOverlap) = 0;
subplot(1,2,2); imagesc(corrMatOverlap);
colormap(redbluecmap);
%%
subplot(2,2,1); imagesc(C); axis square; colorbar
subplot(2,2,2); imagesc(R); axis square; colorbar
subplot(2,2,3); imagesc(R.*C); axis square; colorbar
subplot(2,2,4); imagesc(R.*C>minOverlap); axis square; colorbar
