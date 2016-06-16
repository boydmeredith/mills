peppers = rgb2gray(imread('peppers.png'));
peppers = peppers(1:300,1:300);
yStart = 1;
xStart = 1;
bWidth  = 100;
bHeight = 100;
block   = peppers(yStart:yStart+bHeight-1, xStart:xStart+bWidth-1);
corrMat = normxcorr2(block,peppers);

[yPeak, xPeak] = find(corrMat == max(corrMat(:)));

nbrhdInf.yCtr = yPeak + (bHeight - 1) /2;
nbrhdInf.xMargin = 10;
nbrhdInf.yMargin = 10;

yUL = yPeak - size(block,1) + 1;
xUL = xPeak - size(block,2) + 1;

[corrMatHeight, corrMatWidth] = size(corrMat);

%xCtrInRef = 1-bWidth:20.5:size(peppers,2)+bWidth;
xCtrInRef = linspace(1-(bWidth-1)/2,  size(peppers,2)+(bWidth-1)/2 ,30);
disp(xCtrInRef')
xShiftCorrMat = zeros(length(xCtrInRef), corrMatWidth);


for xx = 1:length(xCtrInRef)
    nbrhdInf.xCtr = xCtrInRef(xx);
    thisCorrMat = computeBlockImageCorrs(block, peppers, nbrhdInf, .5, 'double');
    xShiftCorrMat(xx, :) = thisCorrMat(150,:);
    %xShiftCorrMat(xx, round(xCtrInRef(xx)+bWidth)) = -.5;
end
%
%yCtrInRef = 1-bHeight:20.5:size(peppers,1)+bHeight;
yCtrInRef = linspace(1-(bHeight-1)/2,  size(peppers,1)+(bHeight-1)/2 ,30);

yShiftCorrMat = zeros(corrMatHeight, length(yCtrInRef));

nbrhdInf.xCtr = xPeak + (bWidth - 1) /2;

for yy = 1:length(yCtrInRef)
    nbrhdInf.yCtr = yCtrInRef(yy);
    thisCorrMat = computeBlockImageCorrs(block, peppers, nbrhdInf, .5, 'double');
    yShiftCorrMat(:,yy) = thisCorrMat(:,200);
end

nxc2=normxcorr2(block,peppers);
%
%
figure(5); clf;
subplot(1,2,1)
imagesc([xShiftCorrMat; nxc2(150,:)]); hold on
plot([xCtrInRef(1), xCtrInRef(end)] + (bWidth-1)/2 - max(bWidth,nbrhdInf.xMargin),get(gca, 'ylim'),'k-',...
    [xCtrInRef(1), xCtrInRef(end)] + (3*bWidth-1)/2 - max(bWidth,nbrhdInf.xMargin),get(gca, 'ylim'),'k--',...
    [xCtrInRef(1), xCtrInRef(end)] + (3*bWidth-1)/2 + max(bWidth,nbrhdInf.xMargin),get(gca, 'ylim'),'b--',...
    [xCtrInRef(1), xCtrInRef(end)] + (bWidth-1)/2 + max(bWidth,nbrhdInf.xMargin),get(gca, 'ylim'),'b-')
%set(gca, 'YTick',[1:5:length(xCtrInRef)],'YTickLabel', xCtrInRef(1:5:end), ...
%    'XTick',[1:50:corrMatWidth], 'XTickLabel', [1:50:corrMatWidth] - bWidth + 1);
hold on
plot(xCtrInRef+(bWidth-1)/2,1:size(xShiftCorrMat,1),'ko')
ylabel('x shifts')
xlabel('x in reference')
caxis([-.25 .25])
xlim([1-bWidth/2+1/2, corrMatWidth+bWidth/2-1/2])

subplot(1,2,2)
imagesc(yShiftCorrMat); hold on
plot(get(gca,'xlim'), [yCtrInRef(1), yCtrInRef(end)] + (bHeight-1)/2 - max(bHeight,nbrhdInf.yMargin),'k-',...
    get(gca,'xlim'), [yCtrInRef(1), yCtrInRef(end)] + (3*bHeight-1)/2 - max(bHeight,nbrhdInf.yMargin),'k--',...
    get(gca,'xlim'), [yCtrInRef(1), yCtrInRef(end)] + (3*bHeight-1)/2 + max(bHeight,nbrhdInf.yMargin),'b--',...
    get(gca,'xlim'), [yCtrInRef(1), yCtrInRef(end)] + (bHeight-1)/2 + max(bHeight,nbrhdInf.yMargin),'b-')

hold on
plot(1:size(yShiftCorrMat,2),yCtrInRef+(bHeight-1)/2,'ko')

set(gca,  'YTick',[1:50:corrMatHeight], 'YTickLabel', [1:50:corrMatHeight] - bHeight + 1,...
    'XTick',[1:5:length(yCtrInRef)],'XTickLabel', yCtrInRef(1:5:end))
xlabel('y shifts') 
ylabel('y in reference')
caxis([-.25 .25])
colormap(redbluecmap)


%%




% match  = peppers(yUL:yUL+size(block,1)-1, xUL:xUL+size(block,2)-1);
% mean(block(:)-match(:))
% subplot(1,3,1); imagesc(block);
% subplot(1,3,2); imagesc(match);
% subplot(1,3,3); imagesc(block - match);


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
