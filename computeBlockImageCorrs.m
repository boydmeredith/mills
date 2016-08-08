function [corrMat] = computeBlockImageCorrs(block, reference, nbrhdInf, corrType)
% corrMat = computeBlockImageCorrs(block, reference, corrType); compute the
% normalized 2d cross correlation between a block and a reference. The
% block is expected to be smaller than the reference. 
% 
% inputs:
% block      -      An image selected from a movie frame that we want to
%                   register to the reference
% reference  -      An image in which to look for a match of the block
% nbrhdInf   -      A struct containing information about where we expect
%                   to find the block in the reference with fields
%                       xCtr, yCtr, xMargin, yMargin, minOverlap
% corrType   -      A variable type for storage of the correlation matrix


% if overlap is supplied as a number greater than 1, treat it as a pixel
% count otherwise treat it as a fraction of the block image

[blockHeight, blockWidth] = size(block);
[refHeight, refWidth] = size(reference);
corrMat = zeros(blockHeight+refHeight-1, blockWidth+refWidth-1);
[corrMatHeight, corrMatWidth] = size(corrMat);

if isfield(nbrhdInf,'minOverlap') && nbrhdInf.minOverlap <= 1
    minOverlap = blockHeight * blockWidth * nbrhdInf.minOverlap ; 
elseif isfield(nbrhdInf,'minOverlap')
    minOverlap = ceil(nbrhdInf.minOverlap);
else 
    minOverlap = 0;
end
C = repmat(min(min(1:corrMatWidth,  corrMatWidth:-1:1),blockWidth), corrMatHeight,1);
R = repmat(min(min(1:corrMatHeight, corrMatHeight:-1:1),blockHeight)', 1,corrMatWidth);
legalOverlaps = C.*R >= minOverlap;

if ~isfield(nbrhdInf,'xMargin') || ~isfield(nbrhdInf, 'yMargin') || ...
    (isempty(nbrhdInf.xMargin) && isempty(nbrhdInf.yMargin)) || ...
        ~isfield(nbrhdInf,'xCtr') || ~isfield(nbrhdInf, 'yCtr')
    corrMat = normxcorr2(block, reference);
else
    
    switch 2
        case 1 % enforce a minimum added xmargin of the block dimension + 1
            xMargin  = ceil(max(nbrhdInf.xMargin, blockWidth));
            yMargin  = ceil(max(nbrhdInf.yMargin, blockHeight));
            
        case 2 % no enforcement
            xMargin = nbrhdInf.xMargin;
            yMargin = nbrhdInf.yMargin;
    end
    
    % enforce a block center within the reference image
    %xCtr = min(refWidth, max(1, nbrhdInf.xCtr));
    %yCtr = min(refHeight, max(1, nbrhdInf.yCtr));
    
    % translate from reference image space to correlation matrix space
    xCtrCorr =  floor(nbrhdInf.xCtr + blockWidth/2-1/2);
    yCtrCorr =  floor(nbrhdInf.yCtr + blockHeight/2-1/2);
    
    % put into acceptable zone
    if ~legalOverlaps(round(yCtrCorr),round(xCtrCorr)),
        [yCtrCorr,xCtrCorr] = findClosestAcceptablePoint([yCtrCorr xCtrCorr],legalOverlaps);
    end
    
    % enforce a block center within the safe zone
    %xCtr = min(refWidth, max(1, nbrhdInf.xCtr));
    %yCtr = min(refHeight, max(1, nbrhdInf.yCtr));
    
    % get the offsets of the correlation window that we want to keep
    k1X = max(0, xCtrCorr - 1 - xMargin);
    k1Y = max(0, yCtrCorr - 1 - yMargin);
    k2X = max(0, corrMatWidth - (xCtrCorr + xMargin) - 1);
    k2Y = max(0, corrMatHeight - (yCtrCorr + yMargin)- 1); 
    

    refNbrhd = reference(1+skipNRef(k1Y, blockHeight, refHeight):end-skipNRef(k2Y, blockHeight, refHeight),...
        1+skipNRef(k1X, blockWidth, refWidth):end-skipNRef(k2X, blockWidth, refWidth));
    
    
    %refNbrhd = reference(1+skipNRef(min(k1Y, refHeight - blockHeight + 1), blockHeight):end-skipNRef(min(k2Y, refHeight - blockHeight + 1), blockHeight),...
    %    1+skipNRef(min(k1X, refWidth - blockWidth + 1), blockWidth):end-skipNRef(min(k2X, refWidth - blockWidth + 1), blockWidth));
    
    corrMatNbrhd = normxcorr2(block, refNbrhd);
    
    corrMat(k1Y+1:end-k2Y, k1X+1:end-k2X) = corrMatNbrhd(1+skipNCorr(k1Y, blockHeight, refHeight):end-skipNCorr(k2Y, blockHeight, refHeight),...
        1+skipNCorr(k1X, blockWidth, refWidth):end-skipNCorr(k2X, blockWidth, refWidth));
end

% zero out elements of the correlation matrix that correspond to overlap
% below the minimum required
corrMat = corrMat .* legalOverlaps;

if ~strcmp(corrType,'double')
    corrMat = cast(double(intmax(corrType)) * corrMat,corrType);
end

end


function skipN = skipNRef(k1, B, R)
    skipN = max(min(R-B,1+k1-B),0);
end


function skipN = skipNCorr(k, B,R )
    if 1+k >= B
        if R >= 1+k
            skipN = B-1;
        else
            skipN = k-R+B;
        end
    else
        skipN = k;
    end
end


