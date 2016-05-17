function [block bigBlockRot] = rotateAndReselectBlock(movieFrame, blockLoc, angle);
% block = rotateAndReselectBlock(movieFrame, blockLoc, angle) 
% 
% Rotate movieFrame about the center of the blockLoc by the specified angle
% and return the new image in the blockLoc
%
% movieFrame -          An image from which to select the block 
% blockLoc   -          A logical matrix of the same size as movieFrame with
%                       1's indicating the position of the bloc
% angle      -          An angle in degrees by which to rotate the image
%
% Note: assumes that margins have been built into the edges of the
% blockLocations so that the rotation will not introduce empty pixels to
% the block


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
bigBlockRot = imrotate(movieFrame(bbYPad,bbXPad), angle, 'bilinear' );

[rotH rotW] = size(bigBlockRot);
newDiffH = ceil((rotH - blockH)/2);
newDiffW = ceil((rotW - blockW)/2);
block = bigBlockRot(1+newDiffH : newDiffH+blockH,...
    1+newDiffW : newDiffW+blockW);

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
     

%
%     figure(7); clf
%     imagesc(corrMat);
%     colormap(colormapRedBlue);
end