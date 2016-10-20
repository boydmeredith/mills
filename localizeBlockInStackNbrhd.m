% movieFrame
% thisBlockLoc
% 
% 
% 
% nbrhdInf (xCtr, yCtr, zCtr, rCtr, xMargin, yMargin, rOffToKeep, zOffToKeep)

% params (angleSigFig, rotAngleFromInd, nXYToKeep, corrType)

function [xyzrcoPeak, blockCorrs] = localizeBlockInStackNbrhd(thisBlockLoc,...
    movieFrame, stack, nbrhdInf, varargin)

p = inputParser;
addParameter(p, 'rotAngleFromInd', []);
addParameter(p, 'angleSigFig', 2);
addParameter(p, 'fineRotStepSz', .1);
addParameter(p, 'ignoreZeroRot'  , true);
addParameter(p, 'nXYToKeep'  , Inf);
addParameter(p, 'flagRNFromEdge', 1);
addParameter(p, 'flagZNFromEdge', 1);
addParameter(p, 'whichSlices', []);
p.KeepUnmatched = true;
parse(p,varargin{:});

if ~isfield(nbrhdInf,'fineRotStepSz')
    nbrhdInf.fineRotStepSz = p.Results.fineRotStepSz;
end


if ~isempty(fields(p.Unmatched)),
    error('unmatched fields in localizeBlockInStackNbrhd');
end


whichSlices = p.Results.whichSlices;
if isempty(whichSlices)
    whichSlices = 1:size(stack,3);
end
    

xyzrcoPeak = [];

bInf       = getBlockInf(thisBlockLoc);


% move center of Z search range to closest available stack slice
stackSlices = 1:size(stack,3);
findInVec = @(x,vector)find(abs(x-vector) == min(abs(x-vector)),1,'first');
zCtr = stackSlices(findInVec(nbrhdInf.zCtr,stackSlices));

% choose Z search range to be radius around this center point
zRangeToKeep = intersect(zCtr + nbrhdInf.zOffToKeep,stackSlices);
%zRangeToKeep = zCtr + nbrhdInf.zOffToKeep;
%zRangeToKeep(zRangeToKeep<1 | zRangeToKeep>size(stack,3)) = [];

% choose rotation range to search
rotToKeep = nbrhdInf.rOffToKeep + round(nbrhdInf.rCtr,p.Results.angleSigFig);
rotToKeep = round(rotToKeep, p.Results.angleSigFig);

if ~isempty(p.Results.rotAngleFromInd)
    rotToKeep(~ismember(rotToKeep,p.Results.rotAngleFromInd)) = [];
    blockCorrs.rotAngleFromInd = p.Results.rotAngleFromInd;
else
    blockCorrs.rotAngleFromInd = rotToKeep;
end

if p.Results.ignoreZeroRot
    rotToKeep(rotToKeep==0) = [];
end

if ~isempty(p.Results.nXYToKeep)
    nXYToKeep = 2*(nbrhdInf.xMargin+1) * (nbrhdInf.yMargin + 1);
else
    nXYToKeep = p.Results.nXYToKeep;
end

% create indices for x y z r to keep
nInd = length(zRangeToKeep) * length(rotToKeep) * nXYToKeep;

blockCorrs.indZ = zeros(nInd,1,'uint8');
blockCorrs.indR = zeros(nInd,1,'uint8');
blockCorrs.indX = zeros(nInd,1,'uint16');
blockCorrs.indY = zeros(nInd,1,'uint16');
blockCorrs.corrValsToSave = zeros(nInd,1);
blockCorrs.thisBlockLoc = thisBlockLoc;
blockCorrs.corrMatSize = [size(stack,1)+bInf.height-1 size(stack,2)+bInf.width-1];

thisZRangeToCheck = zRangeToKeep;

nn = 1;
zEdgeFlag = true;

% keep looking for matches until there is a match that is a safe distance
% from the edge of z or we have hit the edge of the stack
while zEdgeFlag && ~isempty(thisZRangeToCheck)
    % iterate through z slices
    for zz = 1:length(thisZRangeToCheck)
        thisSliceNo = thisZRangeToCheck(zz);
        stackSlice = stack(:,:,thisSliceNo);
        % iterate through rotation angles
        for rr = rotToKeep
            % rotate the block according to rr
            blockRot = rotateAndSelectBlock(movieFrame, bInf, rr);
            
            % use find to get y and x indices as well as the values
            % themselves
            [yIx, xIx, thisCorr] = find(computeBlockImageCorrs(blockRot, ...
                stackSlice, nbrhdInf, 'double'));
            
            % don't try to store correlations if we don't find any
            %(bc we are using a uint and all correlations are <= 0)
            if isempty(thisCorr), continue; end
            % --------------------------------------------- %
            
            % get the indices of the top nXYToKeep correlations
            [~, thisCorrSortIx] = sort(thisCorr,'descend');
            
            thisCorrXYToKeepIx  = thisCorrSortIx(1:min(length(thisCorr), nXYToKeep));
            
            % determine where to store these newly computed values
            storeInd = nn:nn+length(thisCorrXYToKeepIx)-1;
            
            % store the x,y,z,r indices and the values themselves for the
            % top nXYToKeep correlations
            blockCorrs.indX(storeInd) = uint16(xIx(thisCorrXYToKeepIx));
            blockCorrs.indY(storeInd) = uint16(yIx(thisCorrXYToKeepIx));
            blockCorrs.indZ(storeInd) = uint8(thisSliceNo);
            blockCorrs.indR(storeInd) = uint8(find(blockCorrs.rotAngleFromInd == rr));
            % store correlations as double (will convert later if desired)
            blockCorrs.corrValsToSave(storeInd) = thisCorr(thisCorrXYToKeepIx);
            
            % increment index counter
            nn = storeInd(end)+1;
        end
    end
    
    
    % get the overall peak for this block/frame
    [cPeak, peakInd] = findPeakCorrVal(blockCorrs, nbrhdInf);
    % if we don't find any peak, get outta here
    if isempty(cPeak), break; end
    
    % check to see if the best z match is too close to an edge
    zTooSmall = find(zRangeToKeep==blockCorrs.indZ(peakInd))-1 <= p.Results.flagZNFromEdge;
    zTooBig  = length(zRangeToKeep)-find(zRangeToKeep==blockCorrs.indZ(peakInd)) <= p.Results.flagZNFromEdge;
    thisZRangeToCheckBig = [];
    thisZRangeToCheckSmall = [];
    if zTooBig,
        thisZRangeToCheckBig = unique(zRangeToKeep(end) + (1:p.Results.flagZNFromEdge));
    end
    if zTooSmall,
        thisZRangeToCheckSmall = unique(zRangeToKeep(1) - (1:p.Results.flagZNFromEdge));
    end
    thisZRangeToCheck = sort( [thisZRangeToCheckBig thisZRangeToCheckSmall]);
    zEdgeFlag = zTooSmall | zTooBig;
    thisZRangeToCheck(~ismember(thisZRangeToCheck,whichSlices)) = [];
    zRangeToKeep = sort([zRangeToKeep thisZRangeToCheck]);
end

% transform x,y from corrMat space to reference space
blkCtrInRefX  = double(blockCorrs.indX(peakInd)) - (bInf.width-1)/2;
blkCtrInRefY  = double(blockCorrs.indY(peakInd)) - (bInf.height-1)/2;
% flag best z and r values that are near edge of computed window
rEdgeFlag = blockCorrs.indR(peakInd)-find(blockCorrs.rotAngleFromInd==rotToKeep(1)) <= p.Results.flagRNFromEdge | ...
    find(blockCorrs.rotAngleFromInd==rotToKeep(end))-blockCorrs.indR(peakInd) <= p.Results.flagRNFromEdge;

if ~isempty(cPeak)
    xyzrcoPeak = [blkCtrInRefX, blkCtrInRefY, double(blockCorrs.indZ(peakInd)), ...
        blockCorrs.rotAngleFromInd(blockCorrs.indR(peakInd)),...
        cPeak, zEdgeFlag | rEdgeFlag]';
end
    
    
end






function [cPeak, peakInd] = findPeakCorrVal(blockCorrs, nbrhdInf)

peakInd = find(blockCorrs.corrValsToSave == max(blockCorrs.corrValsToSave));
cPeak   = max(blockCorrs.corrValsToSave) ;
nPeaks  = length(peakInd);

% if there is only one peak, return it without any further computation
if nPeaks==1, return; end
% if there are more than 10 peaks or the max value is 0, just give up
if nPeaks > 10 || cPeak == 0 || isnan(cPeak), cPeak = []; peakInd = []; return; end


%if nPeaks==1, return, end

pCorrValMean = zeros(nPeaks,1);

for pp = 1:nPeaks
    thisPeakInd = peakInd(pp);
    
    zPNbrhd = blockCorrs.indZ > blockCorrs.indZ(thisPeakInd) - 3 & blockCorrs.indZ < blockCorrs.indZ(thisPeakInd) + 3;
    xPNbrhd = blockCorrs.indX > blockCorrs.indX(thisPeakInd) - 5 & blockCorrs.indX < blockCorrs.indX(thisPeakInd) + 5;
    yPNbrhd = blockCorrs.indY > blockCorrs.indY(thisPeakInd) - 5 & blockCorrs.indY < blockCorrs.indY(thisPeakInd) + 5;
    rPNbrhd = blockCorrs.indR > blockCorrs.indR(thisPeakInd) - 1/nbrhdInf.fineRotStepSz & ...
        blockCorrs.indR < blockCorrs.indR(thisPeakInd) + 1/nbrhdInf.fineRotStepSz;
    
    peakNbrhd = xPNbrhd & yPNbrhd & zPNbrhd & rPNbrhd;disp(sum(peakNbrhd))
    
    pCorrValMean(pp) = mean(blockCorrs.corrValsToSave(peakNbrhd));
    
end

peakInd = peakInd(find(pCorrValMean == max(pCorrValMean)));
% break ties randomly
if length(peakInd) > 1
    warning('WARNING: broke tie between multiple peaks randomly')
    [~, randIx] = max(randn(1,length(peakInd)));
    peakInd = peakInd(randIx);
end
cPeak   = blockCorrs.corrValsToSave(peakInd);
end
