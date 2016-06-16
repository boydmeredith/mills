function [corrValsToSave, xyzrcLocation] = blockwiseMovieStackCorr(subj, movieName, stackName, varargin)

% ---------- parse optional inputs ---------- %
%
p = inputParser;

addOptional(p,'corrType', 'uint16', @(x) ismember(x,{'uint8','uint16','uint32','uint64','double'}));
addOptional(p,'nBlockSpan',10,@(x) isnumeric(x) & ~mod(x,1));
addOptional(p,'blockOverlap',.2,@(x) ispositive(x) & isnumeric(x));

addOptional(p,'coarseRotStepSz',.5, @(x) ispositive(x) & isnumeric(x));
addOptional(p,'coarseRotWindowRange',20, @(x) ispositive(x) & isnumeric(x));
addOptional(p,'fineRotStepSz',.1, @(x) ispositive(x) & isnumeric(x));
addOptional(p,'fineRotWindowRange',2, @(x) ispositive(x) & isnumeric(x));
addOptional(p,'nXYToKeep', 1000, @(x) isnumeric(x) & ~mod(x,1));
addOptional(p,'nZToKeep', 11, @(x) isnumeric(x) & mod(x,2) == 1);
addOptional(p,'zFitPower',5,@(x) isnumeric(x) & ~mod(x,1));
addOptional(p,'angleSigFig',1,@(x) isnumeric(x) & ~mod(x,1));

addOptional(p,'rFitPower',4,@(x) isnumeric & ~mod(x,1));
addOptional(p,'inferZWindow',100,@(x) isnumeric & ~mod(x,1));


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

addOptional(p,'dataDir','/Volumes/tank/jlgauthi/Data',@isdir);

addOptional(p,'summarySaveName', 'summary.mat',@isstr);
addOptional(p,'blockSaveFormat', 'block%03i.mat',@isstr);

addOptional(p, 'rFitFigName','rFit.pdf');
addOptional(p, 'zFitFigName' ,'zFit.pdf');

addOptional(p, 'ballStickGifName', 'ballStick.gif');
addOptional(p, 'diffGifName', 'blockDiffs.gif');
addOptional(p, 'montageGifName', 'blockMontage.gif');
addOptional(p, 'allValsPeakGifName', 'allValsPeak.gif');

addOptional(p, 'priorFigName','priorFig.pdf');

addOptional(p, 'zPriorUseFit', false, @islogical)
addOptional(p, 'rPriorUseFit', true, @islogical)

addOptional(p,'showFigs','off',@(x) any(strcmp({'on','off'},x)));

addOptional(p, 'xyzrPriorSaveName','xyzrPrior.mat',@isstr)
addOptional(p, 'useSavedPrior', true, @islogical);
addOptional(p, 'useSavedPriorEitherWay', false, @islogical);
addOptional(p, 'nbrhdXMargin', 100, @isnumeric);
addOptional(p, 'nbrhdYMargin', 100, @isnumeric);
addOptional(p, 'minCorrOverlap', .5, @isnumeric);
addOptional(p, 'filtKern', [4 4], @(x) isnumeric(x) & length(x)==2)
addOptional(p, 'nRSTD', 8)

%addOptional(p,'saveName','',@isstr);

p.KeepUnmatched = true;

parse(p,varargin{:})
% ----------------------------------------- %

% reassign p.Results to pResults and get rid of p so it is easier to manage
% subfields
pRes = p.Results;
clear p;

%parpool;

if isempty(pRes.nbrhdXMargin) || isempty(pRes.nbrhdYMargin)
    nbrhdInf = [];
else
    nbrhdInf = struct('xMargin', pRes.nbrhdXMargin, 'yMargin', pRes.nbrhdYMargin);
end

set(0, 'DefaultFigureVisible', pRes.showFigs);

close all;

rFitGif  = []; rFitMap  = []; rFitFigNo = 1;
zFitGif  = []; zFitMap  = []; 
ballStickGif   = []; ballStickMap   = []; 
diffGif  = []; diffMap  = []; diffFigNo = 4;
montageGif  = []; montageMap  = []; montageFigNo = 5;
allValsPeakGif = []; allValsPeakMap = [];

if ~isempty(pRes.ballStickGifName)
    ballStickFig = figure('Visible',pRes.showFigs);
end
if ~isempty(pRes.diffGifName)
    diffFig = figure('Visible',pRes.showFigs);
end
if ~isempty(pRes.montageGifName)
    montageFig = figure('Visible',pRes.showFigs);
end
if ~isempty(pRes.allValsPeakGifName)
    allValsPeakFig = figure('Visible',pRes.showFigs);
end

analysisDate = datestr(today);

% assemble some variables based on optional input parameters
pRes.coarseRotAngles  = -(pRes.coarseRotWindowRange/2):pRes.coarseRotStepSz:(pRes.coarseRotWindowRange/2);
pRes.coarseRotAngles = round(pRes.coarseRotAngles, pRes.angleSigFig);
% remove instances of coarseRotAngle == 0 because this will have higher
% correlation than the other angles and screw up our fit
pRes.coarseRotAngles(pRes.coarseRotAngles == 0) = [];
fineRotAngles        = -(pRes.fineRotWindowRange/2):pRes.fineRotStepSz:(pRes.fineRotWindowRange/2);
fineRotAngles        = round(fineRotAngles, pRes.angleSigFig);
rotAngleFromInd      = pRes.coarseRotAngles(1):pRes.fineRotStepSz:pRes.coarseRotAngles(end);
rotAngleFromInd      = round(rotAngleFromInd, pRes.angleSigFig);
pRes.rotAngleFromInd = rotAngleFromInd;
assert(length(pRes.rotAngleFromInd) <= 255); % make sure that we can represent rotation with uint8

nBlocks = pRes.nBlockSpan^2;

% load stack (expect gif)
stackPath = fullfile(pRes.dataDir, subj, stackName);
moviePath = fullfile(pRes.dataDir, subj, movieName);
assert(exist(moviePath,'file') & exist(stackPath,'file'));

stackInf  = imfinfo(stackPath);
stackDim.depth = length(stackInf);

if ~isempty(pRes.loadedStack)
    stack = pRes.loadedStack;
else
    stack = zeros(stackInf(1).Height,stackInf(1).Width,stackDim.depth);
    for ss = 1:length(imfinfo(stackPath))
        stack(:,:,ss) = imread(stackPath,ss);
    end
end
rmfield(pRes,'loadedStack');

% load movie (expect gif)

movieInf = imfinfo(moviePath);
movieHeight = movieInf(1).Height;
movieWidth  = movieInf(1).Width;
movieLength = length(movieInf); 

if ~isempty(pRes.loadedMovie)
    movie = pRes.loadedMovie;
else
    movie = zeros(movieHeight,movieWidth,movieLength);
    for mm = 1:length(imfinfo(moviePath))
        movie(:,:,mm) = imread(moviePath,mm);
    end
end
rmfield(pRes,'loadedMovie');

% check stack and movie size matches file dimensions
assert(isequal([stackInf(1).Height, stackInf(1).Width, stackDim.depth], size(stack)));
assert(isequal([movieHeight, movieWidth, movieLength], size(movie)));

% crop dark edge off stack by removing rows of all zero entries and get new
% size
stack = cropStack(stack);
[stackDim.height, stackDim.width, stackDim.depth] = size(stack);

% create correlations directory
assert(moviePath(end-3)=='.');
movieDate = movieName(1:10);
movieDateDir = fullfile(pRes.dataDir, subj, movieDate);
pRes.corrDir = fullfile(movieDateDir, 'referenceLocalization');

if ~exist(movieDateDir, 'dir'), mkdir(movieDateDir); end
if ~exist(pRes.corrDir, 'dir'), mkdir(pRes.corrDir); end

if isempty(pRes.whichBlocks), pRes.whichBlocks = 1:nBlocks; end
if isempty(pRes.whichFrames), pRes.whichFrames = 1:movieLength; end
if isempty(pRes.whichSlices), pRes.whichSlices = 1:stackDim.depth; end


% divide movie into blocks based on the movie dimensions, desired number of
% blocks, percent overlap between blocks, and maximum amount of rotation
% (used to create a margin)
blockLocations = makeBlockLocations(movieHeight, movieWidth, ...
    pRes.nBlockSpan, pRes.blockOverlap, max(pRes.coarseRotAngles));
maxBWidth = max(max(sum(blockLocations,2),[],1));
maxBHeight = max(max(sum(blockLocations,1),[],2));
montageGrid = zeros(( maxBHeight + 3) * pRes.nBlockSpan, ...
    pRes.nBlockSpan * (3+ maxBWidth));

% initialize matrix to store peak of correlations
xyzrcPeak = zeros(5,nBlocks,movieLength);
xyzrPrior = zeros(5,nBlocks);


% ------ iterate through the frames of the movie ------ %
for ff = 1:length(pRes.whichFrames),
    
    % select the relevant movie frame and the relevant block
    thisFrameNo = pRes.whichFrames(ff);
    movieFrame = movie(:,:,thisFrameNo);
    frameString = sprintf('frame: %03i/%03i',thisFrameNo,movieLength);
    disp(frameString);
    
    % create a directory to store outputs for this frame
    pRes.frameCorrDir = fullfile(pRes.corrDir, sprintf('frame%03i', thisFrameNo));
    if ~exist(pRes.frameCorrDir, 'dir')
        mkdir(pRes.frameCorrDir);
    end
    
    clf(diffFig); clf(allValsPeakFig)
    [~, diffPlotMat]  = makeSubplots(diffFig,sqrt(nBlocks),sqrt(nBlocks),.1,.1,[ 0 0 1 1]);
    [~, allValsPeakMat]  = makeSubplots(allValsPeakFig,sqrt(nBlocks),sqrt(nBlocks),.1,.1,[ 0 0 1 1]);
    
    
    if ff == 1
        xyzrPrior = getPrior(movieFrame, blockLocations, stack, pRes);
    else
        [xx, yy] = meshgrid(1:pRes.nBlockSpan,1:pRes.nBlockSpan);
        [fX, ~, outX]  = fit([xx(:), yy(:)], xyzrcPeak(1,:,thisFrameNo-1)', 'poly11','Robust','Bisquare');
        [fY, ~, outY]  = fit([xx(:), yy(:)], xyzrcPeak(2,:,thisFrameNo-1)', 'poly11','Robust','Bisquare');
        outliersX = outX.residuals > pRes.nRSTD*robustSTD(outX.residuals);
        outliersY = outY.residuals > pRes.nRSTD*robustSTD(outY.residuals);
        outliers = outliersX | outliersY;
        [fZ, ~, outZ] = fit([xx(~outliers), yy(~outliers)], xyzrcPeak(3,~outliers,thisFrameNo-1)', 'loess','Robust','off');
        [fR, ~, outR] = fit([xx(~outliers), yy(~outliers)], xyzrcPeak(4,~outliers,thisFrameNo-1)', 'loess','Robust','off');

        xyzrPrior(1,:) = fX(xx(:),yy(:));
        xyzrPrior(2,:) = fY(xx(:),yy(:));
        xyzrPrior(3,:) = round(fZ(xx(:),yy(:)));
        xyzrPrior(4,:) = round(fR(xx(:),yy(:)),pRes.angleSigFig);
        
    end
    
    % ----------------- iterate through the blocks to keep values ------------------- %
    for bb = 1:length(pRes.whichBlocks)
        % select the relevant block location matrix and get its dimensions so
        % we can determine how many values will be in the correlation matrix
        thisBlockNo   = pRes.whichBlocks(bb);
        thisBlockLoc  = blockLocations(:,:,thisBlockNo);
        bInf          = getBlockInf(thisBlockLoc);
        % get the subscript position for this block in the grid
        [blocki, blockj] = ind2sub([sqrt(nBlocks), sqrt(nBlocks)], thisBlockNo);
        
        if ~isempty(nbrhdInf)
            nbrhdInf.xCtr = xyzrPrior(1,thisBlockNo);
            nbrhdInf.yCtr = xyzrPrior(2,thisBlockNo);
        end
        
        % === Fill in correlations for neighborhood around peak ==== %
        % -----------------------------------------------------------%
        % specify a range of z values and rotation angles to save in the
        % persistent matrix
        thisNbrhdCtrZ = xyzrPrior(3,thisBlockNo);
        thisNbrhdCtrR = xyzrPrior(4,thisBlockNo);
        
        zRangeToKeep = max(1,thisNbrhdCtrZ-ceil(pRes.nZToKeep/2)) : ...
            min(stackDim.depth,thisNbrhdCtrZ+ceil(pRes.nZToKeep/2));
        
        rotToKeep = fineRotAngles + round(thisNbrhdCtrR,pRes.angleSigFig);
        rotToKeep = round(rotToKeep, pRes.angleSigFig);
        rotToKeep(~ismember(rotToKeep,pRes.rotAngleFromInd)) = [];
        
        % create indices for x y z r to keep
        nInd = pRes.nZToKeep * length(rotToKeep) * pRes.nXYToKeep;
        indZ = zeros(nInd,1,'uint8');
        indR = zeros(nInd,1,'uint8');
        indX = zeros(nInd,1,'uint16');
        indY = zeros(nInd,1,'uint16');
        corrValsToSave = zeros(nInd,1,pRes.corrType);
        
        % Loop through the specified range of z values and angles
        fprintf('computing correlations in peak neighborhood for block %03i...\n',thisBlockNo);
        nn = 1;
        for zz = 1:length(zRangeToKeep)
            thisSliceNo = zRangeToKeep(zz);
            stackSlice = stack(:,:,thisSliceNo);
            for rr = rotToKeep
                
                % --- most important lines of the function! --- %
                % rotate the block according to rr
                blockRot = rotateAndSelectBlock(movieFrame, thisBlockLoc, rr);
                
                % use find to get y and x indices as well as the values
                % themselves
                [yIx, xIx, thisCorr] = find(computeBlockImageCorrs(blockRot, ...
                    stackSlice, nbrhdInf, pRes.minCorrOverlap, pRes.corrType));
                % --------------------------------------------- %
                
                % get the indices of the top nXYToKeep correlations
                [~, thisCorrSortIx] = sort(thisCorr,'descend');
                thisCorrXYToKeepIx  = thisCorrSortIx(1:pRes.nXYToKeep);
                
                % determine where to store these newly computed values
                storeInd = nn:nn+pRes.nXYToKeep-1;
                
                % store the x,y,z,r indices and the values themselves for the
                % top nXYToKeep correlations
                indX(storeInd) = uint16(xIx(thisCorrXYToKeepIx));
                indY(storeInd) = uint16(yIx(thisCorrXYToKeepIx));
                indZ(storeInd) = uint8(thisSliceNo);
                indR(storeInd) = uint8(find(pRes.rotAngleFromInd == rr));
                corrValsToSave(storeInd) = thisCorr(thisCorrXYToKeepIx);
                
                % increment index counter
                nn = nn + pRes.nXYToKeep;
            end
        end
        
        % get the peak and save it
        
        [cPeak, peakInd] = findPeakCorrVal(corrValsToSave, indX, indY, indZ, indR, pRes);
        plot(allValsPeakMat(blocki,blockj), 1:length(corrValsToSave), corrValsToSave,'-',...
            peakInd, cPeak, 'mx');
        
        
        blkCtrInRefX  = double(indX(peakInd)) - (bInf.width-1)/2;
        blkCtrInRefY  = double(indY(peakInd)) - (bInf.height-1)/2;
        
        xyzrcPeak(:,thisBlockNo, thisFrameNo) = [blkCtrInRefX, blkCtrInRefY, ...
            double(indZ(peakInd)), pRes.rotAngleFromInd(indR(peakInd)),...
            double(cPeak)/double(intmax(pRes.corrType))]';
        
        
        
        % get the reference block from the peak location. this should
        % probably be its own function
        refBlockIndX = blkCtrInRefX - (bInf.width-1)/2 : blkCtrInRefX + (bInf.width-1)/2;
        refBlockIndY = blkCtrInRefY - (bInf.height-1)/2 : blkCtrInRefY + (bInf.height-1)/2;
        topPad = sum(refBlockIndY < 1);
        bottomPad = sum(refBlockIndY > stackDim.height);
        leftPad = sum(refBlockIndX < 1);
        rightPad = sum(refBlockIndX > stackDim.width);
        
        refBlockIndY(refBlockIndY < 1 | refBlockIndY > stackDim.height) = [];
        refBlockIndX(refBlockIndX < 1 | refBlockIndX > stackDim.width) = [];
        
        refBlock = stack(refBlockIndY, refBlockIndX, indZ(peakInd));
        refBlock = padarray(refBlock, [topPad leftPad], 0, 'pre');
        refBlock = padarray(refBlock, [bottomPad rightPad], 0, 'post');
        
        bestBlockRot = rotateAndSelectBlock(movieFrame, thisBlockLoc, pRes.rotAngleFromInd(indR(peakInd)));
        
        % plot block-ref montage
        montageGrid(1+(maxBHeight+3)*(blocki-1):...
            (maxBHeight+3)*(blocki-1)+size(bestBlockRot,1),...
            1+(2*maxBWidth+3)*(blockj-1):...
            (2*maxBWidth+3)*(blockj-1)+2*size(bestBlockRot,2)) = [normalizeToZeroOne(bestBlockRot) normalizeToZeroOne(refBlock)];
        
        % plot block-ref diff
        imagesc(normalizeToZeroOne(bestBlockRot) - imhistmatch(normalizeToZeroOne(refBlock), bestBlockRot), 'parent',diffPlotMat(blocki,blockj));
        set(diffPlotMat(blocki,blockj), 'XTick',[],'YTick',[]);
        colormap(diffPlotMat(blocki,blockj), bone);
        
        

        
        % % TODO: make two separate functions for pair plots and block
        % znbrhds
        % pair plots
        % %         createPairsPlot(pairsFigNo, peakInd, indX, indY, indZ, indR, ...
        % %             double(corrValsToSave)./double(intmax(pRes.corrType)),...
        % %             frameString, [.2 .9]);
        % %         [pairsGif, pairsMap]  = createGif(pairsFigNo, ff, movieLength, ...
        % %             pairsGif, pairsMap, fullfile(pRes.frameCorrDir, fprintf(pRes.pairsGifName,thisBlockNo)));
        %         % block Z neighbordhood report
        %     [~, blockZNbrhdPlotMat] = makeSubplots(blockZNbrhdFigNo, 4, 7, .1, .1, [ 0 0 1 1]);
        
        %         blockZNbrhdGif = []; blockZNbrhdMap = [];
        %         zNbrhdToPlot = indZ(peakInd)-3:indZ(peakInd)+3;
        %         for zz = 1:length(zNbrhdToPlot)
        %             thisZ = zNbrhdToPlot(zz);
        %             if thisZ < 1 | thisZ > stackDim.depth, continue; end
        %             axes(blockZNbrhdPlotMat(zz,1));
        %             imagesc(refBlock); colormap bone; freezeColors; axis image
        %             ylabel(sprintf('z = %i', thisZ));
        %             axes(blockZNbrhdPlotMat(zz,2));
        %             imagesc(blockRot); colormap bone; freezeColors; axis image
        %             axes(blockZNbrhdPlotMat(zz,3));
        %             imshowpair(refBlock, blockRot); colormap bone; freezeColors; axis image
        %             axes(blockZNbrhdPlotMat(zz,4));
        %
        %         end
        %         [blockZNbrhdGif, blockZNbrhdMap] = createGif(blockZNbrhdFigNo,ff, movieLength,...
        %             blockZNbrhdGif,blockZNbrhdMap,fullfile(pRes.corrDir,blockZNbrhdGifName));
        
        % save block
        if ~isempty(pRes.blockSaveFormat)
            save(fullfile(pRes.frameCorrDir, sprintf(pRes.blockSaveFormat,thisBlockNo)), ...
                'corrValsToSave', 'indX', 'indY', 'indZ', 'indR', 'analysisDate', ...
                'stackPath','rotAngleFromInd','thisBlockLoc');
        end
        
        
        
    end
    
    if ~isempty(pRes.ballStickGifName)
        ballStickFig = xyzrcBallStickFig(xyzrcPeak, thisFrameNo, ballStickFig, stackDim);
        [ballStickGif, ballStickMap] = createGif(ballStickFig,ff,movieLength, ballStickGif,...
            ballStickMap, fullfile(pRes.corrDir,pRes.ballStickGifName));
    end
    
    
    % turn figures into gifs
    %[rFitGif, rFitMap] = createGif(rFitFigNo,ff, movieLength,rFitGif,rFitMap,fullfile(pRes.corrDir,pRes.rFitGifName));
    %[zFitGif, zFitMap] = createGif(zFitFigNo,ff, movieLength,zFitGif,zFitMap,fullfile(pRes.corrDir, pRes.zFitGifName));
    set(0,'currentfigure',montageFig);
    montageFigAx = cla;
    imagesc(montageGrid,'parent', montageFigAx);
    colormap(montageFigAx, 'bone');
    set(montageFig, 'Position', [1, 1, 1577, 954]);
    set(montageFigAx,'position',[.01 .01 .98 .98],'XTick',[],'YTick',[])
    axis(montageFigAx, 'image');
    [montageGif, montageMap] = createGif(montageFig,ff, movieLength,montageGif,...
        montageMap,fullfile(pRes.corrDir, pRes.montageGifName));
    [diffGif, diffMap] = createGif(diffFig,ff, movieLength,diffGif,...
        diffMap,fullfile(pRes.corrDir, pRes.diffGifName));
    [allValsPeakGif, allValsPeakMap] = createGif(allValsPeakFig,ff, movieLength,allValsPeakGif,...
        allValsPeakMap,fullfile(pRes.corrDir, pRes.allValsPeakGifName));
    
    
    % save summary
    if ~isempty(pRes.summarySaveName)
        save(fullfile(pRes.corrDir, pRes.summarySaveName), 'xyzrcPeak', 'blockLocations', ...
            'rotAngleFromInd','stackPath','moviePath','analysisDate');
    end
    
end

set(0, 'DefaultFigureVisible', 'on');

end

%%%%%%% END OF MAIN FUNCTION %%%%%%%%
%========================================================================%
%%%%%%% HELPER FUNCTIONS %%%%%%%%




function [pairsPlot] = createPairsPlot(figNum, peakInd, indX, indY, indZ, indR, corrVals, title, cLimits)
% X-Y
pairsPlot = makeSubplots(figNum,3,2,.1,.1,[0 0 1 1]);
axes(pairsPlot(1));

indZRPeak = find(indZ==indZ(peakInd) & indR ==indR(peakInd));
imagesc(sparse(double(indX(indZRPeak)), double(indY(indZRPeak)), ...
    corrVals(indZRPeak)));
ylabel('Y');
xlabel('X');
axis image; colorbar; colormap(colormapRedBlue); caxis(cLimits);

% X-R
axes(pairsPlot(2));
indZYPeak = find(indZ==indZ(peakInd) & indY ==indY(peakInd));
imagesc(sparse(double(indX(indZYPeak)), double(indR(indZYPeak)), ...
    corrVals(indZYPeak)));
ylabel('R');
xlabel('X');
axis image; colorbar; colormap(colormapRedBlue); caxis(cLimits);

% plot title
text(mean(get(gca,'xlim')), -15,title,'fontsize',25);

% X-Z
axes(pairsPlot(3));
indYRPeak = find(indY==indY(peakInd) & indR ==indR(peakInd));
imagesc(sparse(double(indX(indYRPeak)), double(indZ(indYRPeak)), ...
    corrVals(indYRPeak)));
ylabel('R');
xlabel('X');
axis image; colorbar; colormap(colormapRedBlue); caxis(cLimits);

% Y-R
axes(pairsPlot(4));
indZXPeak = find(indZ==indZ(peakInd) & indX ==indX(peakInd));
imagesc(sparse(double(indY(indZXPeak)), double(indR(indZXPeak)), ...
    corrVals(indZXPeak)));
ylabel('R');
xlabel('Y');
axis image; colorbar; colormap(colormapRedBlue); caxis(cLimits);

% Y-Z
axes(pairsPlot(5));
indXRPeak = find(indX==indX(peakInd) & indR ==indR(peakInd));
imagesc(sparse(double(indY(indXRPeak)), double(indZ(indXRPeak)), ...
    corrVals(indXRPeak)));
xlabel('Y');
ylabel('Z');
axis image; colorbar; colormap(colormapRedBlue); caxis(cLimits);

% Z-R
axes(pairsPlot(6));
indXYPeak = find(indX==indX(peakInd) & indY ==indY(peakInd));
imagesc(sparse(double(indZ(indXYPeak)), double(indR(indXYPeak)), ...
    corrVals(indXYPeak)));
ylabel('R');
xlabel('Z');
axis image; colorbar; colormap(colormapRedBlue); caxis(cLimits);
end

function [cPeak, peakInd] = findPeakCorrVal(corrValsToSave, indX, indY, indZ, indR, pResults)

peakInd = find(corrValsToSave == max(corrValsToSave));
cPeak   = max(corrValsToSave) ;
nPeaks  = length(peakInd);

%if nPeaks==1, return, end

pCorrValMean = zeros(nPeaks,1);

for pp = 1:nPeaks
    thisPeakInd = peakInd(pp);
    
    zPNbrhd = indZ > indZ(thisPeakInd) - 3 & indZ < indZ(thisPeakInd) + 3;
    xPNbrhd = indX > indX(thisPeakInd) - 5 & indX < indX(thisPeakInd) + 5;
    yPNbrhd = indY > indY(thisPeakInd) - 5 & indY < indY(thisPeakInd) + 5;
    rPNbrhd = indR > indR(thisPeakInd) - 1/pResults.fineRotStepSz & ...
        indR < indR(thisPeakInd) + 1/pResults.fineRotStepSz;
    
    peakNbrhd = xPNbrhd & yPNbrhd & zPNbrhd & rPNbrhd;disp(sum(peakNbrhd))
    
    pCorrValMean(pp) = mean(corrValsToSave(peakNbrhd));
    
end

peakInd = peakInd(find(pCorrValMean == max(pCorrValMean)));
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


function [xyzrPrior] = getPrior(movieFrame, blockLocations, stack, pRes)

priorPath = fullfile(pRes.corrDir, pRes.xyzrPriorSaveName);

if pRes.useSavedPrior
    fprintf('Attempting to use saved prior...\n');
    if ~exist(priorPath, 'file')
        fprintf('Could not find file with prior. Recomputing...');
    else
        p = load(priorPath);
        hasSameParams = (p.pRes.nBlockSpan == pRes.nBlockSpan & isequal(p.pRes.whichBlocks, pRes.whichBlocks) ...
            & isequal(p.pRes.whichSlices, pRes.whichSlices) & p.pRes.whichFrames(1) == pRes.whichFrames(1) ...
            & p.pRes.inferZWindow == pRes.inferZWindow & (p.pRes.zFitPower == pRes.zFitPower | pRes.zPriorUseFit == 0) ...
            & p.pRes.zPriorUseFit == pRes.zPriorUseFit &   (p.pRes.rFitPower == pRes.rFitPower | pRes.rPriorUseFit == 0)...
            & p.pRes.rPriorUseFit == pRes.rPriorUseFit & isequal(p.pRes.coarseRotAngles,pRes.coarseRotAngles) ...
            & isequal(p.pRes.rotAngleFromInd,pRes.rotAngleFromInd));
        if hasSameParams || pRes.useSavedPriorEitherWay
            xyzrPrior = p.xyzrPrior;
            return
        else
            fprintf('Parameters did not match previously computed prior. Recomputing...');
            clear p
        end
    end
end

    
if ~isempty(pRes.zFitFigName)
    zFitFig = figure('Visible',pRes.showFigs);
    [~, zFitPlotMat] = makeSubplots(get(zFitFig,'Number'),pRes.nBlockSpan,pRes.nBlockSpan,.1,.1,[.05 .05 .95 .95]);
    linkaxes(zFitPlotMat(:));
end
if ~isempty(pRes.rFitFigName)
    rFitFig = figure('Visible',pRes.showFigs);
    [~, rFitPlotMat] = makeSubplots(get(rFitFig,'Number'),pRes.nBlockSpan,pRes.nBlockSpan,.1,.1,[.05 .05 .95 .95]);
    linkaxes(rFitPlotMat(:));
end

[stackDim.height, stackDim.width, stackDim.depth] = size(stack);
nBlocks = pRes.nBlockSpan^2;

% ------------ get initial z  estimate for each block --------------- %
fprintf('getting initial prior on z for all blocks...\n');
bestZData    = zeros(1,nBlocks);
bestZFit     = zeros(1,nBlocks);
bestRData    = zeros(1,nBlocks);
bestRFit     = zeros(1,nBlocks);
bestXCtrData  = zeros(1,nBlocks);
bestYCtrData  = zeros(1,nBlocks);
blockHeights = zeros(1,nBlocks);
blockWidths  = zeros(1,nBlocks);

for bb = 1:length(pRes.whichBlocks)
    
    % select the relevant block location matrix and get its dimensions so
    % we can determine how many values will be in the correlation matrix
    thisBlockNo  = pRes.whichBlocks(bb);
    thisBlockLoc = blockLocations(:,:,thisBlockNo);
    bInf         = getBlockInf(thisBlockLoc);
    block        = movieFrame(bInf.indY,bInf.indX);
                
    fprintf('computing normxcorr2 for block %03i z estimate...\n',thisBlockNo);
    tic
    
    % matrix to fill to get best z for unrotated blocks
    frameCorrVolNoRot = zeros(bInf.height+stackDim.height-1,bInf.width+stackDim.width-1,...
        stackDim.depth);
    
    for zz = 1:length(pRes.whichSlices),
        thisSliceNo = pRes.whichSlices(zz);
        stackSlice = stack(:,:,thisSliceNo);
        frameCorrVolNoRot(:,:,thisSliceNo) = computeBlockImageCorrs(block, ...
            stackSlice, [], pRes.minCorrOverlap, 'double');
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
    
    zIndToFit   = max(min(pRes.whichSlices), bestZData(thisBlockNo)-ceil(pRes.inferZWindow/2)) : ...
        min(max(pRes.whichSlices), bestZData(thisBlockNo)+ceil(pRes.inferZWindow/2));

    fitZ = polyfit(double(zIndToFit), double(maxCorrByZData(zIndToFit)'), pRes.zFitPower);
    
    maxCorrByZFit = polyval(fitZ, zIndToFit);
    bestZFit(thisBlockNo) = zIndToFit(maxCorrByZFit == max(maxCorrByZFit));
    
    % save a report of the peak correlations and fits that will be used to get neighborhoods 
    if ~isempty(pRes.zFitFigName)
        [blocki, blockj] = ind2sub([pRes.nBlockSpan, pRes.nBlockSpan], thisBlockNo);
        
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
        
        if blockj == 1 & blocki == ceil(pRes.nBlockSpan/2), ylabel(ax1,'correlation');
        end
        if blockj == ceil(pRes.nBlockSpan/2) & blocki == pRes.nBlockSpan, xlabel(ax1,'stack slice #')
        end
    end
   
    clear frameCorrVolNoRot;
end

if ~isempty(pRes.zFitFigName), saveas(zFitFig, fullfile(pRes.frameCorrDir, pRes.zFitFigName)); end

[xx, yy] = meshgrid(1:pRes.nBlockSpan,1:pRes.nBlockSpan);

[fX, ~, outX]  = fit([xx(:), yy(:)], bestXCtrData', 'poly11','Robust','Bisquare');
[fY, ~, outY]  = fit([xx(:), yy(:)], bestYCtrData', 'poly11','Robust','Bisquare');
outliersX = outX.residuals > pRes.nRSTD*robustSTD(outX.residuals);
outliersY = outY.residuals > pRes.nRSTD*robustSTD(outY.residuals);
outliers = outliersX | outliersY;

if pRes.zPriorUseFit
    [fZ, ~, outZ] = fit([xx(~outliers), yy(~outliers)], bestZFit(~outliers)', 'loess','Robust','off');
else
    [fZ, ~, outZ] = fit([xx(~outliers), yy(~outliers)], bestZData(~outliers)', 'loess','Robust','off');
end
xyzrPrior(1,:) = fX(xx(:),yy(:));
xyzrPrior(2,:) = fY(xx(:),yy(:));
xyzrPrior(3,:) = round(fZ(xx(:),yy(:)));
assert(isequal(fZ(xx(:),yy(:)),reshape(fZ(xx,yy),[],1)));

% ------------ get initial rotation angle estimate for each block --------------- %
fprintf('getting initial prior on r for all blocks...\n');
highestCorr=-Inf;lowestCorr=Inf;
for bb = 1:length(pRes.whichBlocks)
    % select the relevant block location matrix and get its dimensions so
    % we can determine how many values will be in the correlation matrix
    thisBlockNo  = pRes.whichBlocks(bb);
    thisBlockLoc = blockLocations(:,:,thisBlockNo);
    bInf         = getBlockInf(thisBlockLoc);
    block        = movieFrame(bInf.indY,bInf.indX);
    
    fprintf('computing normxcorr2 for block %03i r estimate...\n',thisBlockNo);
    tic
    
    % use the block z prior that we just computed above
    thisBlockBestZ = xyzrPrior(3,thisBlockNo);
    bestStackSliceNoRot = stack(:,:,thisBlockBestZ);
    
    % preallocate matrix for rotation angles
    frameCorrVolRot = zeros(bInf.height+stackDim.height-1,bInf.width+stackDim.width-1,...
        length(pRes.coarseRotAngles), 'double');
    
    % get correlations for a pre-determined set of angles
    for rr = 1:length(pRes.coarseRotAngles)
        rotAngle = pRes.coarseRotAngles(rr);
        blockRot = rotateAndSelectBlock(movieFrame, thisBlockLoc, rotAngle);
        frameCorrVolRot(:,:,rr) = computeBlockImageCorrs(blockRot, ...
            bestStackSliceNoRot, [], pRes.minCorrOverlap, 'double');
    end
    toc
    
    % find peak correlation for rotations using same process as for z
    [~, frameCorrVolMaxIx] = max(frameCorrVolRot(:));
    [bestYData, bestXData, bestRIndData] =  ind2sub(size(frameCorrVolRot), frameCorrVolMaxIx);

    % figure out what angle corresponds to the best r index 
    bestRData(thisBlockNo) = pRes.coarseRotAngles(bestRIndData);
    maxCorrByRData = squeeze(max(max(frameCorrVolRot,[],1),[],2));
    
    fitRot =  polyfit((pRes.coarseRotAngles).',maxCorrByRData, pRes.rFitPower);
    maxCorrByRFit = polyval(fitRot,pRes.rotAngleFromInd);
    
    [~, bestRIndFit] = max(maxCorrByRFit);
    bestRFit(thisBlockNo) = pRes.rotAngleFromInd(bestRIndFit);
    
    if highestCorr < max(maxCorrByRData), highestCorr = max(maxCorrByRData); end
    if lowestCorr > min(maxCorrByRData), lowestCorr   = min(maxCorrByRData); end

    % save a report of the peak correlations and fits that will be used to get neighborhoods 
    if ~isempty(pRes.rFitFigName)
        [blocki, blockj] = ind2sub([pRes.nBlockSpan, pRes.nBlockSpan], thisBlockNo);
        ax1 = rFitPlotMat(blocki, blockj);
        hold(ax1, 'on');
        set(ax1,'box','off','ylim',[lowestCorr highestCorr], 'YTick',[round(max(maxCorrByRData),2)],...%'YTickLabel',[max(maxCorrByRData)], ...
            'xlim',[pRes.coarseRotAngles(1) pRes.coarseRotAngles(end)], 'XTick',[bestRData(thisBlockNo)],...
            'fontsize', 6); 
        
        % plot fits
        plot(ax1,repmat(bestRFit(thisBlockNo),1,2),[0 1], '-.' ...
            ,pRes.rotAngleFromInd, maxCorrByRFit, '-', 'color', [.9 .5 .7], 'linewidth',1);
        % plot data 
        plot(ax1,repmat(bestRData(thisBlockNo),1,2),[0 1], '--', 'color', [.3 .5 .7],'linewidth',1)
        plot(ax1, pRes.coarseRotAngles, maxCorrByRData, '.', 'color', [.2 .4 .6],'markersize',1.75);
        
        if blockj == 1 & blocki == ceil(pRes.nBlockSpan/2), ylabel(ax1,'correlation');
        end
        if blockj == ceil(pRes.nBlockSpan/2) & blocki == pRes.nBlockSpan, xlabel(ax1,'rotation angle (deg)')
        end
    end
    
    clear frameCorrVolRot;
    
end

if ~isempty(pRes.rFitFigName), saveas(rFitFig, fullfile(pRes.frameCorrDir, pRes.rFitFigName)); end
[fX, ~, outX]  = fit([xx(:), yy(:)], bestXCtrData', 'poly11','Robust','Bisquare');
[fY, ~, outY]  = fit([xx(:), yy(:)], bestYCtrData', 'poly11','Robust','Bisquare');
outliersX = outX.residuals > pRes.nRSTD*robustSTD(outX.residuals);
outliersY = outY.residuals > pRes.nRSTD*robustSTD(outY.residuals);
outliers = outliersX | outliersY;

if pRes.zPriorUseFit
    [fZ, ~, outZ] = fit([xx(~outliers), yy(~outliers)], bestZFit(~outliers)', 'loess','Robust','off');
else
    [fZ, ~, outZ] = fit([xx(~outliers), yy(~outliers)], bestZData(~outliers)', 'loess','Robust','off');
end
xyzrPrior(1,:) = fX(xx(:),yy(:));
xyzrPrior(2,:) = fY(xx(:),yy(:));
xyzrPrior(3,:) = round(fZ(xx(:),yy(:)));
assert(isequal(fZ(xx(:),yy(:)),reshape(fZ(xx,yy),[],1)));

if pRes.zPriorUseFit
    [fR, ~, outR] = fit([xx(~outliers), yy(~outliers)], bestRFit(~outliers)', 'loess','Robust','off');
else
    [fR, ~, outR] = fit([xx(~outliers), yy(~outliers)], bestRData(~outliers)', 'loess','Robust','off');
end
xyzrPrior(4,:) = round(fR(xx(:),yy(:)),pRes.angleSigFig);
if pRes.zPriorUseFit
    [fR, ~, outR] = fit([xx(~outliers), yy(~outliers)], bestRFit(~outliers)', 'loess','Robust','off');
else
    [fR, ~, outR] = fit([xx(~outliers), yy(~outliers)], bestRData(~outliers)', 'loess','Robust','off');
end
xyzrPrior(4,:) = round(fR(xx(:),yy(:)),pRes.angleSigFig);


if ~isempty(pRes.priorFigName)
    priorFig = figure('visible', pRes.showFigs);
    
    [~, pfM]      = makeSubplots(priorFig, 4, 2, .1, .1, [.05 .05 .95 .95]);
    colormap(priorFig, colormapRedBlue);
    
    imagesc([reshape(bestXCtrData,pRes.nBlockSpan,pRes.nBlockSpan)...
        reshape(xyzrPrior(1,:),pRes.nBlockSpan,pRes.nBlockSpan)],'parent',pfM(1,1));
    title(pfM(1,1),'X');     colorbar(pfM(1,1))
    axis(pfM(1,1),'image');

    imagesc(reshape(outliersX,pRes.nBlockSpan,pRes.nBlockSpan),'parent',pfM(2,1));
    
    hist(pfM(2,2), outX.residuals,1000);
    hold(pfM(2,2), 'on')
    plot(pfM(2,2), pRes.nRSTD*[robustSTD(outX.residuals) robustSTD(outX.residuals)],get(pfM(2,2),'ylim'),'r',...
    -pRes.nRSTD*[robustSTD(outX.residuals) robustSTD(outX.residuals)],get(pfM(2,2),'ylim'),'r')
    
    imagesc([reshape(bestYCtrData,pRes.nBlockSpan,pRes.nBlockSpan)...
        reshape(xyzrPrior(2,:),pRes.nBlockSpan,pRes.nBlockSpan)],'parent',pfM(1,2));
    title(pfM(1,2),'Y'); colorbar(pfM(1,2))
    axis(pfM(1,2),'image');
    
    imagesc(reshape(outliersY,pRes.nBlockSpan,pRes.nBlockSpan),'parent',pfM(2,3));
    
    hist(pfM(2,4), outY.residuals,1000);
    hold(pfM(2,4), 'on')
    plot(pfM(2,4), pRes.nRSTD*[robustSTD(outY.residuals) robustSTD(outY.residuals)],get(pfM(2,4),'ylim'),'r',...
    -pRes.nRSTD*[robustSTD(outY.residuals) robustSTD(outY.residuals)],get(pfM(2,4),'ylim'),'r')

    imagesc([reshape(bestZData,pRes.nBlockSpan,pRes.nBlockSpan)...
        reshape(xyzrPrior(3,:),pRes.nBlockSpan,pRes.nBlockSpan)],'parent',pfM(1,3));
    title(pfM(1,3),'Z'); colorbar(pfM(1,3))
    axis(pfM(1,3),'image');
    
    imagesc([reshape(bestRData,pRes.nBlockSpan,pRes.nBlockSpan)...
        reshape(xyzrPrior(4,:),pRes.nBlockSpan,pRes.nBlockSpan)],'parent',pfM(1,4));
    title(pfM(1,4),'R'); colorbar(pfM(1,4))
    axis(pfM(1,4),'image');
    
    set(priorFig, 'position',[1 1 1600 1000], 'paperpositionmode','auto');
    saveas(priorFig, fullfile(pRes.corrDir, pRes.priorFigName));
    close(priorFig)
end

save(priorPath, 'xyzrPrior','pRes')

end

