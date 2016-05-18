function [blockLocations] = makeBlockLocations(imgH, imgW, nBars, overlap, maxRot)
% [blockLocations] = makeBlockLocations(imgH, imgW, nBars, percentOverlap, maxRot)
% Divide an image into partially overlapping square blocks using Jeff's
% function divideIntoBlocks and return the block locations
% imgH              height of image to divvy up
% imgW              weidth of image to divvy up
% nBars             the number of blocks to span the image vertically and
%                       horizontally
% overlap           amount of overlap to give the blocks 
%                       if >=1, assume # points
%                       otherwise, assume fractional overlap
% maxRot            an amount of rotation that we should allow for by
%                   padding the edges of the image (see intaglio doc for
%                   more info)
%
% blockLocations    a binary matrix for each block showing its location

canRotate  = false;
imgHOrig = imgH;
imgWOrig = imgW;
borderW = 0;
borderH = 0;

while ~canRotate
    
    [startsH, endsH, lengthH] = divideIntoBlocks(imgH, nBars, overlap);
    [startsW, endsW, lengthW] = divideIntoBlocks(imgW, nBars, overlap);

    blockLocations = false(imgH, imgW, nBars^2);
    bb = 1;
    for bj = 1:nBars
        for bi = 1:nBars %
            blockLocations([startsH(bi):endsH(bi)],...
                [startsW(bj):endsW(bj)],bb) = true;

            bb = bb+1;

        end
    end

    blockW = endsW(1);
    blockH = endsH(1);
    % test ability to rotate within original image
    % if we rotate a block, how much wider and taller will it be?
    [padH padW] = dimAfterRotation(blockH, blockW, maxRot);
    rotDiffW = ceil(padW - blockW);
    rotDiffH = ceil(padH - blockH);
    
    canRotate = ( borderW >= rotDiffW/2 & borderH >= rotDiffH/2 );
   
    borderW = borderW + rotDiffW;
    borderH = borderH + rotDiffH;
    imgH = imgHOrig - borderH;
    imgW = imgWOrig - borderW;
    
    
end

% pad blockLocations with borderW and borderH
blockLocations = padarray(blockLocations,[ceil(borderH/2) ceil(borderW/2)],0);

