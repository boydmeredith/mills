blockWidth = 50;
blockHeight = 60;
overlap = .5;
minOverlap = blockHeight * blockWidth * overlap ; 
if overlap > 1
    minOverlap = overlap;
end
corrMatWidth = 512;
corrMatHeight = 510;

C = repmat(min(1:corrMatWidth,corrMatWidth:-1:1),corrMatHeight,1);
R = repmat(min(1:corrMatHeight, corrMatHeight:-1:1)',1,corrMatWidth);

C = min(C, blockWidth);
R = min(R, blockHeight);
subplot(2,2,1); imagesc(C); axis square; colorbar
subplot(2,2,2); imagesc(R); axis square; colorbar
subplot(2,2,3); imagesc(R.*C); axis square; colorbar
subplot(2,2,4); imagesc(R.*C>minOverlap); axis square; colorbar
