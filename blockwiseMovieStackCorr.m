function [corrValsToSave, xyzrcoPeak] = blockwiseMovieStackCorr(subj, movieDate, varargin)

% ---------- parse optional inputs ---------- %
%
p = inputParser;

addOptional(p, 'location','L01');

addOptional(p, 'stackDate', []);

addOptional(p,'corrType', 'uint16', @(x) ismember(x,{'uint8','uint16','uint32','uint64','double'}));
addOptional(p,'mByNBlocks',[6 6],@(x) isnumeric(x) & ~mod(x,1));
%addOptional(p,'blockOverlap',.2,@(x) ispositive(x) & isnumeric(x));
addOptional(p,'blockOverlap',10,@(x) ispositive(x) & isnumeric(x));

addOptional(p,'coarseRotStepSz',.5, @(x) ispositive(x) & isnumeric(x));
addOptional(p,'coarseRotWindowRange',20, @(x) ispositive(x) & isnumeric(x));

if strcmp(subj,'J114')
    addOptional(p,'fineRotWindowRange',10, @(x) ispositive(x) & isnumeric(x));
    addOptional(p,'fineRotStepSz',.5, @(x) ispositive(x) & isnumeric(x));
    
else
    addOptional(p,'fineRotWindowRange',6, @(x) ispositive(x) & isnumeric(x));
    addOptional(p,'fineRotStepSz',.25, @(x) ispositive(x) & isnumeric(x));
    
end
addOptional(p,'nXYToKeep', 400, @(x) isnumeric(x) & ~mod(x,1));

addOptional(p,'nZToKeepInlier', 11, @(x) isnumeric(x) & mod(x,2) == 1);
addOptional(p,'nZToKeepOutlier', 21, @(x) isnumeric(x) & mod(x,2) == 1);

%addOptional(p,'angleSigFig',2,@(x) isnumeric(x) & ~mod(x,1));
addOptional(p,'angleSigFig',2,@(x) isnumeric(x) & ~mod(x,1));

% inferZWindow is only relevant if useZFit
addOptional(p,'inferZWindow',100,@(x) isnumeric & ~mod(x,1));
addOptional(p, 'zSearchRangeUseFit', false, @islogical)
addOptional(p, 'rSearchRangeUseFit', false, @islogical)
addOptional(p,'zFitPower',5,@(x) isnumeric(x) & ~mod(x,1));
addOptional(p,'rFitPower',4,@(x) isnumeric & ~mod(x,1));


addOptional(p,'whichBlocks', [],@isnumeric);
addOptional(p,'whichFrames', [],@isnumeric);
addOptional(p,'whichSlices', [],@isnumeric);

addOptional(p,'reportZPeakFinder',true,@islogical);
addOptional(p,'reportRotPeakFinder',true,@islogical);

addOptional(p,'reportBlockRefZNbrhdd',true,@islogical);
addOptional(p,'reportBlockRefRotNbrhdd',true,@islogical);

addOptional(p,'reportBlockRefOverlap',true,@islogical);
addOptional(p,'reportBlockRefDiff',true,@islogical);
addOptional(p,'reportDimPairs',true,@islogical);

addOptional(p,'loadedStack',[],@isnumeric);
addOptional(p,'loadedMovie',[],@isnumeric);

addOptional(p,'dataDir',[],@(x) isdir(x) | isempty(x));

addOptional(p,'summarySaveName', 'summary.mat',@isstr);
addOptional(p,'blockSaveFormat', 'block%03i.mat',@isstr);

addOptional(p, 'rFitFigName','rFit.pdf');
addOptional(p, 'zFitFigName' ,'zFit.pdf');

addOptional(p, 'diffGifName', 'blockDiffs.gif');
addOptional(p, 'montageGifName', 'blockMontage.gif');
addOptional(p, 'allValsPeakGifName', 'allValsPeak.gif');

addOptional(p, 'searchRangeFigName','searchRangeFig.pdf');



addOptional(p,'showFigs','off',@(x) any(strcmp({'on','off'},x)));

addOptional(p, 'xyzrSearchRangeSaveName','xyzrSearchRange.mat',@isstr)
addOptional(p, 'useSavedSearchRange', true, @islogical);
addOptional(p, 'useSavedSearchRangeEitherWay', false, @islogical);
addOptional(p, 'nbrhdXMargin', 10, @isnumeric);
addOptional(p, 'nbrhdYMargin', 10, @isnumeric);
addOptional(p, 'searchRangeXMargin', 10, @isnumeric);
addOptional(p, 'searchRangeYMargin', 10, @isnumeric);
addOptional(p, 'minCorrOverlap', .8, @isnumeric);
addOptional(p, 'nRSTD', 8)
addOptional(p,'xRadiusMin',4,@isnumeric);
addOptional(p,'yRadiusMin',4,@isnumeric);

addOptional(p, 'flagZNFromEdge', 2)
addOptional(p, 'flagRNFromEdge', 2)



%addOptional(p,'saveName','',@isstr);

p.KeepUnmatched = true;

parse(p,varargin{:})
% ----------------------------------------- %


fprintf('\nstarting registration for subj: %s on movie from %s \n',subj, movieDate);
% reassign p.Results to params and get rid of p so it is easier to manage
% subfields
if ~isempty(fields(p.Unmatched))
    warning(sprintf('unmatched fields in input: \n\t%s\n', strjoin(fields(p.Unmatched),'\n\t')));
end
params = p.Results;
clear p;

%parpool;

params.subj = subj;
params.movieDate = movieDate;
if isempty(params.dataDir)
    params.dataDir = jlgDataDir;
end



if isempty(params.nbrhdXMargin) || isempty(params.nbrhdYMargin)
    nbrhdInf = [];
else
    nbrhdInf = struct('xMargin', params.nbrhdXMargin, 'yMargin', params.nbrhdYMargin);
end

% if only one number is given for mByNBlocks, repeat it for both dimensions
if size(params.mByNBlocks) == 1, params.mByNBlocks = repmat(params.mByNBlocks, 1, 2); end
params.nBlocks = prod(params.mByNBlocks);

if (params.nbrhdYMargin+1)*2 * (params.nbrhdXMargin+1)*2 < params.nXYToKeep,
    params.nXYToKeep = (params.nbrhdYMargin+1)*2 * (params.nbrhdXMargin+1)*2;
end

set(0, 'DefaultFigureVisible', params.showFigs);

% date information to save
dateStr = datestr(now);
dateNum = now;

% assemble some variables based on optional input parameters
params.coarseRotAngles  = -(params.coarseRotWindowRange/2):params.coarseRotStepSz:(params.coarseRotWindowRange/2);
params.coarseRotAngles = round(params.coarseRotAngles, params.angleSigFig);
% remove instances of coarseRotAngle == 0 because this will have higher
% correlation than the other angles and screw up our fit
params.coarseRotAngles(params.coarseRotAngles == 0) = [];
fineRotAngles          = -(params.fineRotWindowRange/2): params.fineRotStepSz ...
    :(params.fineRotWindowRange/2);
fineRotAngles          = round(fineRotAngles, params.angleSigFig);
rotAngleFromInd        = params.coarseRotAngles(1):params.fineRotStepSz ...
    :params.coarseRotAngles(end);
rotAngleFromInd        = round(rotAngleFromInd, params.angleSigFig);
params.rotAngleFromInd = rotAngleFromInd;
assert(length(params.rotAngleFromInd) <= 255); % make sure that we can represent rotation with uint8



% load stack (expect gif)
movieFname = sprintf('%s__%s__AVERAGE.tif',movieDate,params.location);
if isempty(params.stackDate),
    params.stackDate =  defaultStackDate(subj);
end
stackPath = fullfile(subj, sprintf('reference_stack_%s.tif',params.stackDate));
fullStackPath = fullfile(params.dataDir, stackPath);
moviePath = fullfile(subj, movieFname);
fullMoviePath = fullfile(params.dataDir, moviePath);
assert(exist(fullMoviePath,'file') & exist(fullStackPath,'file'));

stackInf  = imfinfo(fullStackPath);
stackDim.depth = length(stackInf);

if ~isempty(params.loadedStack)
    stack = params.loadedStack;
else
    stack = zeros(stackInf(1).Height,stackInf(1).Width,stackDim.depth);
    for ss = 1:length(imfinfo(fullStackPath))
        stack(:,:,ss) = imread(fullStackPath,ss);
    end
end
rmfield(params,'loadedStack');

% load movie (expect gif)

movieInf = imfinfo(fullMoviePath);
movieHeight = movieInf(1).Height;
movieWidth  = movieInf(1).Width;
movieLength = length(movieInf);

if ~isempty(params.loadedMovie)
    movie = params.loadedMovie;
else
    movie = zeros(movieHeight,movieWidth,movieLength);
    for mm = 1:length(imfinfo(fullMoviePath))
        movie(:,:,mm) = imread(fullMoviePath,mm);
    end
end
rmfield(params,'loadedMovie');

% check stack and movie size matches file dimensions
assert(isequal([stackInf(1).Height, stackInf(1).Width, stackDim.depth], size(stack)));
assert(isequal([movieHeight, movieWidth, movieLength], size(movie)));

% crop dark edge off stack by removing rows of all zero entries and get new
% size
stack = cropStack(stack);
[stackDim.height, stackDim.width, stackDim.depth] = size(stack);

% create correlations directory
assert(fullMoviePath(end-3)=='.');
movieDate = movieFname(1:10);
movieDateDir = fullfile(params.dataDir, subj, movieDate);
params.corrDir = referenceLocalizationDir(subj, movieDate, params.location);

if ~exist(movieDateDir, 'dir'), mkdir(movieDateDir); end
if ~exist(params.corrDir, 'dir'), mkdir(params.corrDir); end

if isempty(params.whichBlocks), params.whichBlocks = 1:params.nBlocks; end
if isempty(params.whichFrames), params.whichFrames = 1:movieLength; end
if isempty(params.whichSlices), params.whichSlices = 1:stackDim.depth; end


% divide movie into blocks based on the movie dimensions, desired number of
% blocks, percent overlap between blocks, and maximum amount of rotation
% (used to create a margin)
blockLocations = makeBlockLocations(movieHeight, movieWidth, ...
    params.mByNBlocks, params.blockOverlap, max(params.coarseRotAngles));


% initialize matrix to store peak of correlations
xyzrcoPeak = zeros(6,params.nBlocks,movieLength);
xyzrSearchRange = zeros(5,params.nBlocks);

if ~isempty(params.summarySaveName)
    summaryPath = fullfile(params.corrDir, params.summarySaveName);
    save(summaryPath,'xyzrcoPeak', 'blockLocations','rotAngleFromInd',...
    'stackPath','moviePath','dateStr', 'dateNum','params',...
    '-append');
    summfile = matfile(summaryPath,'writable',true);
end

% ------ iterate through the frames of the movie ------ %
for ff = 1:length(params.whichFrames),
    
    % select the relevant movie frame and the relevant block
    thisFrameNo = params.whichFrames(ff);
    movieFrame = movie(:,:,thisFrameNo);
    frameString = sprintf('frame: %03i/%03i',thisFrameNo,movieLength);
    disp(frameString);
    
    % create a directory to store outputs for this frame
    params.frameCorrDir = fullfile(params.corrDir, sprintf('frame%03i', thisFrameNo));
    if ~exist(params.frameCorrDir, 'dir')
        mkdir(params.frameCorrDir);
    end
    
    if ff == 1
        [xyzrSearchRange, outliersXY] = getSearchRange(movieFrame, blockLocations, stack, params);
    end
    
    % ----------------- iterate through the blocks to keep values ------------------- %
    for bb = 1:length(params.whichBlocks)
        % select the relevant block location matrix and get its dimensions so
        % we can determine how many values will be in the correlation matrix
        thisBlockNo   = params.whichBlocks(bb);
        thisBlockLoc  = blockLocations(:,:,thisBlockNo);
        bInf          = getBlockInf(thisBlockLoc);
        % compute size of full correlation matrix based on block info
        corrMatSize = [bInf.height + stackDim.height - 1, bInf.width + stackDim.width - 1];
        
        % === Fill in correlations for neighborhood around peak ==== %
        % -----------------------------------------------------------%
        if ~isempty(nbrhdInf)
            nbrhdInf.xCtr = xyzrSearchRange(1,thisBlockNo);
            nbrhdInf.yCtr = xyzrSearchRange(2,thisBlockNo);
        end
        
        % specify a range of z values and rotation angles to save in the
        % persistent matrix
        
        thisNbrhdCtrZ = xyzrSearchRange(3,thisBlockNo);
        thisNbrhdCtrR = xyzrSearchRange(4,thisBlockNo);
        
        if outliersXY(thisBlockNo),
            thisNZToKeep = params.nZToKeepOutlier;
        else
            thisNZToKeep = params.nZToKeepInlier;
        end
        
        zRangeToKeep = max(1,thisNbrhdCtrZ-ceil(thisNZToKeep/2)) : ...
            min(stackDim.depth,thisNbrhdCtrZ+ceil(thisNZToKeep/2));
        thisZRangeToCheck = zRangeToKeep;
        
        rotToKeep = fineRotAngles + round(thisNbrhdCtrR,params.angleSigFig);
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
        
        % Loop through the specified range of z values and angles
        fprintf('computing correlations in peak neighborhood for block %03i...\n',thisBlockNo);
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
        else
            xyzrcoPeak(:, thisBlockNo, thisFrameNo) = [xyzrSearchRange(:,thisBlockNo); nan; 1];
        end
        
        if ~isempty(params.blockSaveFormat)
            save(fullfile(params.frameCorrDir, sprintf(params.blockSaveFormat,thisBlockNo)), ...
                'corrValsToSave', 'indX', 'indY', 'indZ', 'indR', 'dateStr','dateNum', ...
                'stackPath','rotAngleFromInd','thisBlockLoc','corrMatSize');
        end
    end
    
    % update the search range for the next frame and keep track of the
    % outliers for this frame
    [xyzrSearchRange, outliersXY] = fitXYZRSearchRange(xyzrcoPeak(1,:, thisFrameNo),xyzrcoPeak(2,:, thisFrameNo),...
        xyzrcoPeak(3,:, thisFrameNo), xyzrcoPeak(4,:, thisFrameNo), params);
    xyzrcoPeak(6,:, thisFrameNo) = xyzrcoPeak(6,:, thisFrameNo) | outliersXY;
    
    % write the results of this frame into the summary file
    if ~isempty(params.summarySaveName)
        summfile.xyzrcoPeak(:,:,thisFrameNo) = xyzrcoPeak(:,:,thisFrameNo);
    end
    
end

set(0, 'DefaultFigureVisible', 'on');

end







%%%%%%% END OF MAIN FUNCTION %%%%%%%%
%========================================================================%
%%%%%%% HELPER FUNCTIONS %%%%%%%%


function [cPeak, peakInd] = findPeakCorrVal(corrValsToSave, indX, indY, indZ, indR, params)

peakInd = find(corrValsToSave == max(corrValsToSave));
cPeak   = max(corrValsToSave) ;
nPeaks  = length(peakInd);

% if there is only one peak, return it without any further computation
if nPeaks==1, return; end
% if there are more than 10 peaks or the max value is 0, just give up
if nPeaks > 10 || cPeak == 0 || isnan(cPeak), cPeak = []; peakInd = []; return; end


%if nPeaks==1, return, end

pCorrValMean = zeros(nPeaks,1);

for pp = 1:nPeaks
    thisPeakInd = peakInd(pp);
    
    zPNbrhd = indZ > indZ(thisPeakInd) - 3 & indZ < indZ(thisPeakInd) + 3;
    xPNbrhd = indX > indX(thisPeakInd) - 5 & indX < indX(thisPeakInd) + 5;
    yPNbrhd = indY > indY(thisPeakInd) - 5 & indY < indY(thisPeakInd) + 5;
    rPNbrhd = indR > indR(thisPeakInd) - 1/params.fineRotStepSz & ...
        indR < indR(thisPeakInd) + 1/params.fineRotStepSz;
    
    peakNbrhd = xPNbrhd & yPNbrhd & zPNbrhd & rPNbrhd;disp(sum(peakNbrhd))
    
    pCorrValMean(pp) = mean(corrValsToSave(peakNbrhd));
    
end

peakInd = peakInd(find(pCorrValMean == max(pCorrValMean)));
% break ties randomly
if length(peakInd) > 1
    warning('WARNING: broke tie between multiple peaks randomly')
    [~, randIx] = max(randn(1,length(peakInd)));
    peakInd = peakInd(randIx);
end
cPeak   = corrValsToSave(peakInd);
end

function [peaks yPeakInRef xPeakInRef rPeakInRef] = getPeakCorrCoords(corrMat, blksz, rotationAngles)
peaks = ind2sub(size(corrMat),find(corrMat==max(corrMat(:))));
yPeakInRef = yPeak - blksz + 1;
xPeakInRef = xPeak - blksz + 1;
rPeakInRef = rotationAngles(rPeak);
end


function [starts ends] = getRangeToCheck(bestIx, maxIx, nToKeep)
% determine the window of different z values to look in
starts  = bestIx-(nToKeep-1)/2;
ends    = bestIx+(nToKeep-1)/2;
if starts < 1
    ends = ends-(starts-1);
    starts = 1;
elseif ends > maxIx
    starts = starts-(ends-maxIx);
    ends = maxIx;
end
end


function [xyzrSearchRange, outliersXY] = getSearchRange(movieFrame, blockLocations, stack, params)
% function [xyzrSearchRange] = getSearchRange(movieFrame, blockLocations, stack, params)

if isempty(params.searchRangeXMargin) || isempty(params.searchRangeYMargin)
    nbrhdInf = [];
else
    nbrhdInf = struct('xMargin', params.searchRangeXMargin, 'yMargin', params.searchRangeYMargin);
end

if ~isempty(params.summarySaveName)
    summaryPath = fullfile(params.corrDir, params.summarySaveName);
    % load the summary.mat file
    sfile = matfile(summaryPath, 'Writable',true);
    if ~isfield(sfile,'params')
        fprintf('Could not find xyzrSearchRange in summary.mat. Recomputing...');
    else
        fprintf('Attempting to use saved searchRange...\n');
        hasSameParams = (isequal(sfile.params.mByNBlocks, params.mByNBlocks) & isequal(sfile.params.whichBlocks, params.whichBlocks) ...
            & isequal(sfile.params.whichSlices, params.whichSlices) & sfile.params.whichFrames(1) == params.whichFrames(1) ...
            & sfile.params.inferZWindow == params.inferZWindow & (sfile.params.zFitPower == params.zFitPower | params.zSearchRangeUseFit == 0) ...
            & sfile.params.zSearchRangeUseFit == params.zSearchRangeUseFit &   (sfile.params.rFitPower == params.rFitPower | params.rSearchRangeUseFit == 0)...
            & sfile.params.rSearchRangeUseFit == params.rSearchRangeUseFit & isequal(sfile.params.coarseRotAngles,params.coarseRotAngles) ...
            & isequal(sfile.params.rotAngleFromInd,params.rotAngleFromInd) & params.nRSTD==sfile.params.nRSTD & params.xRadiusMin==sfile.params.xRadiusMin & ...
            params.yRadiusMin==sfile.params.yRadiusMin);

        if (hasSameParams || params.useSavedSearchRangeEitherWay) && isfield(sfile,'xyzrSearchRange') && isfield(sfile,'outliersXY')
            xyzrSearchRange = sfile.xyzrSearchRange;
            outliersXY      = sfile.outliersXY;
            return
        else
            fprintf('Parameters did not match previously computed searchRange. Recomputing...');
        end
    end
end


if ~isempty(params.zFitFigName)
    zFitFig = figure('Visible',params.showFigs);
    [~, zFitPlotMat] = makeSubplots(get(zFitFig,'Number'),params.mByNBlocks(2),params.mByNBlocks(1),.1,.1,[.05 .05 .95 .95]);
    linkaxes(zFitPlotMat(:));
end
if ~isempty(params.rFitFigName)
    rFitFig = figure('Visible',params.showFigs);
    [~, rFitPlotMat] = makeSubplots(get(rFitFig,'Number'),params.mByNBlocks(2),params.mByNBlocks(1),.1,.1,[.05 .05 .95 .95]);
    linkaxes(rFitPlotMat(:));
end

[stackDim.height, stackDim.width, stackDim.depth] = size(stack);


% ------------ get initial z  estimate for each block --------------- %
fprintf('getting initial searchRange on z for all blocks...\n');
bestZData    = zeros(1,params.nBlocks);
bestZFit     = zeros(1,params.nBlocks);
bestRData    = zeros(1,params.nBlocks);
bestRFit     = zeros(1,params.nBlocks);
bestXCtrData  = zeros(1,params.nBlocks);
bestYCtrData  = zeros(1,params.nBlocks);
bestXCtrDataRot = zeros(1,params.nBlocks);
bestYCtrDataRot = zeros(1,params.nBlocks);
blockHeights = zeros(1,params.nBlocks);
blockWidths  = zeros(1,params.nBlocks);

for bb = 1:length(params.whichBlocks)
    
    % select the relevant block location matrix and get its dimensions so
    % we can determine how many values will be in the correlation matrix
    thisBlockNo  = params.whichBlocks(bb);
    thisBlockLoc = blockLocations(:,:,thisBlockNo);
    bInf         = getBlockInf(thisBlockLoc);
    block        = movieFrame(bInf.indY,bInf.indX);
    
    fprintf('computing normxcorr2 for block %03i z estimate...\n',thisBlockNo);
    tic
    
    % matrix to fill to get best z for unrotated blocks
    frameCorrVolNoRot = zeros(bInf.height+stackDim.height-1,bInf.width+stackDim.width-1,...
        stackDim.depth);
    
    for zz = 1:length(params.whichSlices),
        thisSliceNo = params.whichSlices(zz);
        stackSlice = stack(:,:,thisSliceNo);
        frameCorrVolNoRot(:,:,thisSliceNo) = computeBlockImageCorrs(block, ...
            stackSlice, [], params.minCorrOverlap, 'double');
    end
    toc
    
    % find maximum correlation for each z and fit a polynomial to the
    % correlation values in a window around the max. Then choose z that
    % maximizes the fit curve
    
    maxCorrByZData = squeeze(max(max(frameCorrVolNoRot(:,:,:),[],1),[],2));
    % get position of upper left corner of block at best match
    [~, maxInd]    = max(frameCorrVolNoRot(:));
    [bestYData, bestXData, bestZData(thisBlockNo)] = ind2sub(size(frameCorrVolNoRot), maxInd);
    
    bestYCtrData(thisBlockNo) = bestYData - bInf.height/2 + 1/2;
    bestXCtrData(thisBlockNo) = bestXData - bInf.width/2 + 1/2;
    blockHeights(thisBlockNo) = bInf.height;
    blockWidths(thisBlockNo)  = bInf.width;
    
    zIndToFit   = max(min(params.whichSlices), bestZData(thisBlockNo)-ceil(params.inferZWindow/2)) : ...
        min(max(params.whichSlices), bestZData(thisBlockNo)+ceil(params.inferZWindow/2));
    
    fitZ = polyfit(double(zIndToFit), double(maxCorrByZData(zIndToFit)'), params.zFitPower);
    
    maxCorrByZFit = polyval(fitZ, zIndToFit);
    bestZFit(thisBlockNo) = zIndToFit(maxCorrByZFit == max(maxCorrByZFit));
    
    % save a report of the peak correlations and fits that will be used to get neighborhoods
    if ~isempty(params.zFitFigName)
        [blocki, blockj] = ind2sub(params.mByNBlocks, thisBlockNo);
        
        ax1 = zFitPlotMat(blocki, blockj); hold(ax1, 'on');
        set(ax1,'box','off','ylim',[-.25 1], 'YTick',[round(max(maxCorrByZData),2)],...
            'xlim',[1 stackDim.depth], 'XTick', bestZData(thisBlockNo) , 'fontsize',6);
        
        % plot fit to data
        plot(ax1,[bestZFit(thisBlockNo) bestZFit(thisBlockNo)],[0 1], '-.', ...
            zIndToFit, maxCorrByZFit, '-', 'color', [.9 .5 .7],'linewidth',1);
        % plot data points with mark best z
        plot(ax1,repmat(bestZData(thisBlockNo),1,2),[0 1], '--','color', [.3 .5 .7],'linewidth',1);
        plot(ax1, 1:stackDim.depth, maxCorrByZData, '.', 'color', [.2 .4 .6], 'markersize',1.75)%;,...,
        %'markeredgecolor','k','linewidth',.05,'markersize',3);
        
        if blockj == 1 & blocki == ceil(params.mByNBlocks(1)/2), ylabel(ax1,'correlation');
        end
        if blockj == ceil(params.mByNBlocks(2)/2) & blocki == params.mByNBlocks(1), xlabel(ax1,'stack slice #')
        end
    end
    
    clear frameCorrVolNoRot;
end

if ~isempty(params.zFitFigName), saveas(zFitFig, fullfile(params.frameCorrDir, params.zFitFigName)); end

% fit x, y, z to get search range from correlations obtained so far
if params.zSearchRangeUseFit
    zIn = bestZFit;
else
    zIn = bestZData;
end
[xyzSearchRange, outliersXY] = fitXYZRSearchRange(bestXCtrData,bestYCtrData,zIn,[],params);

% ------------ get initial rotation angle estimate for each block --------------- %
fprintf('getting initial searchRange on r for all blocks...\n');
highestCorr=-Inf;lowestCorr=Inf;
for bb = 1:length(params.whichBlocks)
    % select the relevant block location matrix and get its dimensions so
    % we can determine how many values will be in the correlation matrix
    thisBlockNo  = params.whichBlocks(bb);
    thisBlockLoc = blockLocations(:,:,thisBlockNo);
    bInf         = getBlockInf(thisBlockLoc);
    block        = movieFrame(bInf.indY,bInf.indX);
    
    fprintf('computing normxcorr2 for block %03i r estimate...\n',thisBlockNo);
    tic
    
    % use the block z searchRange that we just computed above
    bestStackSliceNoRot = stack(:,:,xyzSearchRange(3,thisBlockNo));
    % set up x and y ctr
    if ~isempty(nbrhdInf)
        nbrhdInf.xCtr = xyzSearchRange(1,thisBlockNo);
        nbrhdInf.yCtr = xyzSearchRange(2,thisBlockNo);
    end
    
    
    
    % preallocate matrix for rotation angles
    frameCorrVolRot = zeros(bInf.height+stackDim.height-1,bInf.width+stackDim.width-1,...
        length(params.coarseRotAngles), 'double');
    
    % get correlations for a pre-determined set of angles
    for rr = 1:length(params.coarseRotAngles)
        rotAngle = params.coarseRotAngles(rr);
        
        blockRot = rotateAndSelectBlock(movieFrame, bInf, rotAngle);
        frameCorrVolRot(:,:,rr) = computeBlockImageCorrs(blockRot, ...
            bestStackSliceNoRot, nbrhdInf, params.minCorrOverlap, 'double');
        
    end
    toc
    
    % find peak correlation for rotations using same process as for z
    [~, frameCorrVolMaxIx] = max(frameCorrVolRot(:));
    [bestYDataRot, bestXDataRot, bestRIndData] =  ind2sub(size(frameCorrVolRot), frameCorrVolMaxIx);
    bestXCtrDataRot(thisBlockNo) = bestXDataRot - bInf.width/2 + 1/2;
    bestYCtrDataRot(thisBlockNo) = bestYDataRot - bInf.height/2 + 1/2;
    % figure out what angle corresponds to the best r index
    bestRData(thisBlockNo) = params.coarseRotAngles(bestRIndData);
    maxCorrByRData = squeeze(max(max(frameCorrVolRot,[],1),[],2));
    
    fitRot =  polyfit((params.coarseRotAngles).',maxCorrByRData, params.rFitPower);
    maxCorrByRFit = polyval(fitRot,params.rotAngleFromInd);
    
    [~, bestRIndFit] = max(maxCorrByRFit);
    bestRFit(thisBlockNo) = params.rotAngleFromInd(bestRIndFit);
    
    % keep track of best and worst correlations in order to set axis limits
    % later
    if highestCorr < max(maxCorrByRData), highestCorr = max(maxCorrByRData); end
    if lowestCorr > min(maxCorrByRData), lowestCorr   = min(maxCorrByRData); end
    
    % save a report of the peak correlations and fits that will be used to get neighborhoods
    if ~isempty(params.rFitFigName)
        [blocki, blockj] = ind2sub(params.mByNBlocks, thisBlockNo);
        ax1 = rFitPlotMat(blocki, blockj);
        hold(ax1, 'on');
        set(ax1,'box','off','ylim',[lowestCorr highestCorr], 'YTick',[round(max(maxCorrByRData),2)],...%'YTickLabel',[max(maxCorrByRData)], ...
            'xlim',[params.coarseRotAngles(1) params.coarseRotAngles(end)], 'XTick',[bestRData(thisBlockNo)],...
            'fontsize', 6);
        
        % plot fits
        plot(ax1,repmat(bestRFit(thisBlockNo),1,2),[0 1], '-.' ...
            ,params.rotAngleFromInd, maxCorrByRFit, '-', 'color', [.9 .5 .7], 'linewidth',1);
        % plot data
        plot(ax1,repmat(bestRData(thisBlockNo),1,2),[0 1], '--', 'color', [.3 .5 .7],'linewidth',1)
        plot(ax1, params.coarseRotAngles, maxCorrByRData, '.', 'color', [.2 .4 .6],'markersize',1.75);
        
        if blockj == 1 & blocki == ceil(params.mByNBlocks(1)/2), ylabel(ax1,'correlation');
        end
        if blockj == ceil(params.mByNBlocks(2)/2) & blocki == params.mByNBlocks(1), xlabel(ax1,'rotation angle (deg)')
        end
    end
    
    clear frameCorrVolRot;
    
end

if ~isempty(params.rFitFigName), saveas(rFitFig, fullfile(params.frameCorrDir, params.rFitFigName)); end



% fit x, y, z, r to get search range from correlations obtained looking at
% rotation
if params.rSearchRangeUseFit
    rIn = bestRFit;
else
    rIn = bestRData;
end
[xyzrSearchRange, outliersXY] = fitXYZRSearchRange(bestXCtrDataRot,bestYCtrDataRot,...
    xyzSearchRange(3,:),rIn,params,true);

if ~isempty(params.summarySaveName)
    % write to the summary.mat file
    sfile.xyzrSearchRange             = xyzrSearchRange;
    sfile.xyzrSearchRangeOutliersXY   = outliersXY;
end

end






