function [blockLocations] = makeBlockLocations(img, nBars, overlap)
% [blockLocations] = getBlockLocations(image, nBars, percentOverlap)
% Divide an image into partially overlapping square blocks using Jeff's
% function divideIntoBlocks and return the block locations
% image             the image to divvy up
% nBars             the number of blocks to span the image vertically and
%                       horizontally
% overlap           amount of overlap to give the blocks 
%                       if >=1, assume # points
%                       otherwise, assume fractional overlap
%
% blockLocations    a binary matrix for each block showing its location

[imgH, imgW] = size(img); 
[startsH, endsH, lengthH] = divideIntoBlocks(imgH, nBars, overlap);
[startsW, endsW, lengthW] = divideIntoBlocks(imgW, nBars, overlap);

blockLocations = zeros(imgH, imgW, nBars^2);
bb = 1;
for wb = 1:nBars
    for hb = 1:nBars
        blockLocations([startsH(hb):endsH(hb)],...
            [startsW(wb):endsW(wb)],bb) = 1;
        bb = bb+1;
        
    end
end