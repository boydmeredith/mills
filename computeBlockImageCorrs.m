function [corrMat block] = computeBlockImageCorrs(movieFrame, blockLoc, angle, reference, corrType);
% corrMat = computeBlockImageCorrs(movieFrame, reference, blockLoc, angle,
% corrType); compute the normalized 2d cross correlation between a block
% from the image movieFrame (with location specified by the logical matrix
% blockLoc) after rotating the block by the specified angle.
% movieFrame -      An image from which to select the block
% blockLoc   -      A logical matrix of the same size as movieFrame with
%                   1's indicating the position of the bloc
% angle      -      An angle in degrees by which to rotate the image
% reference  -      An image in which to look for a match of the block
% corrType   -      A variable type for storage of the correlation matrix

% get size of the movie and of the block
[movieH, movieW] = size(movieFrame);
[bbY, bbX] = getBlockIx(blockLoc);
blockW = length(bbX);
blockH = length(bbY);

% determine how large a window around the block we need in order to be able
% to rotate the block around its center by the specified angle without
% introducing empty pixels into the block
[padH padW] = dimAfterRotation(blockH, blockW, angle);

halfRotDiffW = ceil((padW - blockW)/2);
halfRotDiffH = ceil((padH - blockH)/2);

bbYPad = (bbY(1) - halfRotDiffH) : (bbY(end) + halfRotDiffH) ;
bbYPad = bbYPad(bbYPad > 0 & bbYPad <= movieH);
bbXPad = (bbX(1) - halfRotDiffW) : (bbX(end) + halfRotDiffW);
bbXPad = bbXPad(bbXPad > 0 & bbXPad <= movieW);

% rotate the block and reselect that chunk that we really want
bigBlockRot = imrotate(movieFrame(bbYPad,bbXPad), angle);

[rotH rotW] = size(bigBlockRot);
newDiffH = ceil((rotH - blockH)/2);
newDiffW = ceil((rotW - blockW)/2);
block = bigBlockRot(1+newDiffH : newDiffH+blockH,...
    1+newDiffW : newDiffH+blockW);

%     figure(11); clf
%     subplot(2,2,1)
%     imagesc(movieFrame(bbY,bbX)); axis image
%     subplot(2,2,2)
%     imagesc(movieFrame(bbYPad,bbXPad)); axis image
%     subplot(2,2,3)
%     imagesc(bigBlockRot); axis image
%     subplot(2,2,4)
%     imagesc(block); axis image
%     figure(11)
     

corrMat = normxcorr2(block,reference);
if ~strcmp(corrType,'double')
    corrMat = cast(double(intmax(corrType)) * corrMat,corrType);
end
%
%     figure(7); clf
%     imagesc(corrMat);
%     colormap(colormapRedBlue);
end