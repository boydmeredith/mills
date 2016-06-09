function [corrMat] = computeBlockImageCorrs(block, reference, nbrhdInf, minOverlap, corrType)
% corrMat = computeBlockImageCorrs(block, reference, corrType); compute the
% normalized 2d cross correlation between a block and a reference. The
% block is expected to be smaller than the reference. 
% 
% inputs:
% block      -      An image selected from a movie frame that we want to
%                   register to the reference
% reference  -      An image in which to look for a match of the block
% minOverlap -      Specifies the minimal amount of overlap between the
%                   block and the reference to keep in the correlation
%                   matrix
% nbrhdInf   -      A struct containing information about where we expect
%                   to find the block in the reference with fields
%                       xCtr, yCtr, xMargin, yMargin
% corrType   -      A variable type for storage of the correlation matrix


% if overlap is supplied as a number greater than 1, treat it as a pixel
% count otherwise treat it as a fraction of the block image

[blockHeight, blockWidth] = size(block);
[refHeight, refWidth] = size(reference);
corrMat = zeros(blockHeight+refHeight-1, blockWidth+refWidth-1);
[corrMatHeight, corrMatWidth] = size(corrMat);

if minOverlap <= 1
    minOverlap = blockHeight * blockWidth * minOverlap ; 
else
    minOverlap = ceil(minOverlap);
end
C = min(repmat(min(1:corrMatWidth,  corrMatWidth:-1:1),corrMatHeight,1),   blockWidth);
R = min(repmat(min(1:corrMatHeight, corrMatHeight:-1:1)',1,corrMatWidth),blockHeight);
legalOverlaps = C.*R > minOverlap;

if isempty(nbrhdInf)
    corrMat = normxcorr2(block, reference);
else
    % enforce a minimum added xmargin of the block dimension + 1
    xMargin  = ceil(max(nbrhdInf.xMargin, blockWidth));
    yMargin  = ceil(max(nbrhdInf.yMargin, blockHeight));
    % enforce a block center within the reference image
    xCtr = min(refWidth, max(1, nbrhdInf.xCtr));
    yCtr = min(refHeight, max(1, nbrhdInf.yCtr));
    % get the offsets of the correlation window that we want to keep
    xStartCorr =  floor(xCtr + blockWidth/2-1/2);
    yStartCorr =  floor(yCtr + blockHeight/2-1/2);
    k1X = max(0, xStartCorr - 1 - xMargin);
    k1Y = max(0, yStartCorr - 1 - yMargin);
    k2X = max(0, corrMatWidth - (xStartCorr + xMargin) - 1);
    k2Y = max(0, corrMatHeight - (yStartCorr + yMargin)- 1); 
    
    
%     minOverlapWidth  = minOverlap/blockHeight;
%     minOverlapHeight = minOverlap/blockWidth;
%     nbrhdXCtr   = min(max(minOverlapWidth, nbrhdInf.xCtr), refWidth-minOverlapWidth);
%     nbrhdYCtr   = min(max(minOverlapHeight, nbrhdInf.yCtr), refHeight-minOverlapHeight);
    
    refNbrhd = reference(1+skipNRef(k1Y, blockHeight):end-skipNRef(k2Y, blockHeight),...
        1+skipNRef(k1X, blockWidth):end-skipNRef(k2X, blockWidth));
    corrMatNbrhd = normxcorr2(block, refNbrhd);
    
    corrMat(k1Y+1:end-k2Y, k1X+1:end-k2X) = corrMatNbrhd(1+skipNCorr(k1Y, blockHeight):end-skipNCorr(k2Y, blockHeight),...
        1+skipNCorr(k1X, blockWidth):end-skipNCorr(k2X, blockWidth));
end

% zero out elements of the correlation matrix that correspond to overlap
% below the minimum required
corrMat(~legalOverlaps) = 0;

if ~strcmp(corrType,'double')
    corrMat = cast(double(intmax(corrType)) * corrMat,corrType);
end

end


function skipN = skipNRef(k1, B)
    skipN = max(1+k1-B, 0);
end


function skipN = skipNCorr(k, B)
    if 1+k >= B
        skipN = B-1;
    else
        skipN = k;
    end
end


