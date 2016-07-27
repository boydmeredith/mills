function [xyzrcoPeak, thisBlockCorrs] = blockwiseMovieStackCorrFinalStep(bInf, movieFrame, stack, nbrhdInf,params)
if outliersXY(thisBlockNo),
    thisNZToKeep = params.nZToKeepOutlier;
else
    thisNZToKeep = params.nZToKeepInlier;
end
% compute size of full correlation matrix based on block info
corrMatSize = [bInf.height + stackDim.height - 1, bInf.width + stackDim.width - 1];

zRangeToKeep = max(1,nbrhdInf.ctrZ-ceil(thisNZToKeep/2)) : ...
    min(stackDim.depth,nbrhdInf.ctrZ+ceil(thisNZToKeep/2));
thisZRangeToCheck = zRangeToKeep;

rotToKeep = fineRotAngles + round(nbrhdInf.ctrR,params.angleSigFig);
rotToKeep = round(rotToKeep, params.angleSigFig);
rotToKeep(~ismember(rotToKeep,params.rotAngleFromInd)) = [];
rotToKeep(rotToKeep==0) = [];

% create indices for x y z r to keep
nInd = thisNZToKeep * length(rotToKeep) * params.nXYToKeep;
indZ = zeros(nInd,1,'uint8');
indR = zeros(nInd,1,'uint8');
indX = zeros(nInd,1,'uint16');
indY = zeros(nInd,1,'uint16');
corrValsToSave = zeros(nInd,1,params.corrType);

nn = 1;
zEdgeFlag = true;
while zEdgeFlag && ~isempty(thisZRangeToCheck)
    for zz = 1:length(thisZRangeToCheck)
        thisSliceNo = thisZRangeToCheck(zz);
        stackSlice = stack(:,:,thisSliceNo);
        for rr = rotToKeep
            
            % --- most important lines of the function! --- %
            % rotate the block according to rr
            blockRot = rotateAndSelectBlock(movieFrame, bInf, rr);
            
            % use find to get y and x indices as well as the values
            % themselves
            [yIx, xIx, thisCorr] = find(computeBlockImageCorrs(blockRot, ...
                stackSlice, nbrhdInf, params.minCorrOverlap, 'double'));
            
            % don't try to store correlations if we don't find any
            %(bc we are using a uint and all correlations are <= 0)
            if isempty(thisCorr), continue; end
            % --------------------------------------------- %
            
            % get the indices of the top nXYToKeep correlations
            [~, thisCorrSortIx] = sort(thisCorr,'descend');
            
            thisCorrXYToKeepIx  = thisCorrSortIx(1:min(length(thisCorr), params.nXYToKeep));
            
            % determine where to store these newly computed values
            storeInd = nn:nn+length(thisCorrXYToKeepIx)-1;
            
            % store the x,y,z,r indices and the values themselves for the
            % top nXYToKeep correlations
            indX(storeInd) = uint16(xIx(thisCorrXYToKeepIx));
            indY(storeInd) = uint16(yIx(thisCorrXYToKeepIx));
            indZ(storeInd) = uint8(thisSliceNo);
            indR(storeInd) = uint8(find(params.rotAngleFromInd == rr));
            % convert correlations to desired type and store them
            if strcmp(params.corrType,'double')
                corrValsToSave(storeInd) = thisCorr(thisCorrXYToKeepIx);
            else
                corrValsToSave(storeInd) = cast(thisCorr(thisCorrXYToKeepIx) *...
                    double(intmax(params.corrType)), params.corrType);
            end
            
            % increment index counter
            nn = storeInd(end)+1;
        end
    end
    
    
    % get the overall peak for this block/frame
    [cPeak, peakInd] = findPeakCorrVal(corrValsToSave, indX, indY, indZ, indR, params);
    % if we don't find any peak, get outta here
    if isempty(cPeak), break; end
    
    % check to see if the best z match is too close to an edge
    zTooSmall = find(zRangeToKeep==indZ(peakInd))-1 <= params.flagZNFromEdge;
    zTooBig  = length(zRangeToKeep)-find(zRangeToKeep==indZ(peakInd)) <= params.flagZNFromEdge;
    if zTooBig,
        thisZRangeToCheck = zRangeToKeep(end) + (1:params.flagZNFromEdge);
    elseif zTooSmall,
        thisZRangeToCheck = zRangeToKeep(1) - (1:params.flagZNFromEdge);
    end
    zEdgeFlag =  zTooSmall | zTooBig ;
    thisZRangeToCheck(~ismember(thisZRangeToCheck,params.whichSlices)) = [];
    zRangeToKeep = sort([zRangeToKeep thisZRangeToCheck]);
end

% transform x,y from corrMat space to reference space
blkCtrInRefX  = double(indX(peakInd)) - (bInf.width-1)/2;
blkCtrInRefY  = double(indY(peakInd)) - (bInf.height-1)/2;
% flag best z and r values that are near edge of computed window
rEdgeFlag = indR(peakInd)-find(params.rotAngleFromInd==rotToKeep(1)) <= params.flagRNFromEdge | ...
    find(params.rotAngleFromInd==rotToKeep(end))-indR(peakInd) <= params.flagRNFromEdge;

% save absolute peak
if ~isempty(cPeak)
    if ~strcmp(params.corrType,'double'),
        cPeak = double(cPeak)/double(intmax(params.corrType));
    end
    xyzrcoPeak(:,thisBlockNo, thisFrameNo) = [blkCtrInRefX, blkCtrInRefY, ...
        double(indZ(peakInd)), params.rotAngleFromInd(indR(peakInd)),...
        cPeak, zEdgeFlag | rEdgeFlag]';
else % if no good peak found, use the center of the search range
    xyzrcoPeak(:, thisBlockNo, thisFrameNo) = [xyzrSearchRange(:,thisBlockNo); nan; 1];
    % round the x and y centers according to the dimensions of the
    % block. If the block dimension is even, the center should be at a
    % multiple of .5
    thisX = xyzrcoPeak(1, thisBlockNo, thisFrameNo);
    thisY = xyzrcoPeak(2, thisBlockNo, thisFrameNo);
    % if width is even, put x center at closest .5, else closest integer
    xyzrcoPeak(1, thisBlockNo, thisFrameNo) = round(thisX) + ~mod(bInf.width,2)*sign(thisX-round(thisX))*.5;
    % if height is even, put y center at closest .5, else closest integer
    xyzrcoPeak(2, thisBlockNo, thisFrameNo) = round(thisY) + ~mod(bInf.height,2)*sign(thisY-round(thisY))*.5;
end

% save block-specific registration information
if ~isempty(params.blockSaveFormat)
    save(fullfile(params.frameCorrDir, sprintf(params.blockSaveFormat,thisBlockNo)), ...
        'corrValsToSave', 'indX', 'indY', 'indZ', 'indR', 'dateStr','dateNum', ...
        'stackPath','rotAngleFromInd','thisBlockLoc','corrMatSize');
end
end

