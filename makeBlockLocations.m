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
rotDiffW = 0;
rotDiffH = 0;

while ~canRotate

    borderW = borderW + rotDiffW + mod(rotDiffW,2);
    borderH = borderH + rotDiffH + mod(rotDiffH,2);
    
    imgH = imgHOrig - borderH;
    imgW = imgWOrig - borderW;
    
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
    
    % get dimensions of biggest block
    maxBlockW = max(sum(max(blockLocations,[],1)));
    maxBlockH = max(sum(max(blockLocations,[],2)));
    
    % test ability to rotate within original image
    % if we rotate a block, how much wider and taller will it be?
    [padH padW] = dimAfterRotation(maxBlockH, maxBlockW, maxRot);

    rotDiffW = ceil(padW) - maxBlockW;
    rotDiffH = ceil(padH) - maxBlockH;

    canRotate = ( borderW > rotDiffW & borderH > rotDiffH );
        
    
    
end

assert(imgH + borderH == imgHOrig);
assert(imgW + borderW == imgWOrig);

% pad blockLocations with borderW and borderH
blockLocations = padarray(blockLocations,[(borderH/2) (borderW/2)],0);
assert(isequal(size(blockLocations), [imgHOrig imgWOrig nBars^2]));
