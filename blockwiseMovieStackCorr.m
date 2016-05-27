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
addOptional(p,'nXYToKeep', 30, @(x) isnumeric(x) & ~mod(x,1));
addOptional(p,'nZToKeep', 11, @(x) isnumeric(x) & mod(x,2) == 1);
addOptional(p,'zFitPower',5,@(x) isnumeric(x) & ~mod(x,1));
addOptional(p,'angleSigFig',1,@(x) isnumeric(x) & ~mod(x,1));

addOptional(p,'rotFitPower',4,@(x) isnumeric & ~mod(x,1));
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

addOptional(p, 'priorMedFiltFigName','priorMedFilt.pdf');

addOptional(p,'showFigs','off',@(x) any(strcmp({'on','off'},x)));

%addOptional(p,'saveName','',@isstr);

p.KeepUnmatched = true;

parse(p,varargin{:})
% ----------------------------------------- %

% reassign p.Results to pResults and get rid of p so it is easier to manage
% subfields
pRes = p.Results;
clear p;

%parpool;

set(0, 'DefaultFigureVisible', pRes.showFigs);

close all;

rFitGif  = []; rFitMap  = []; rFitFigNo = 1;
zFitGif  = []; zFitMap  = []; 
ballStickGif   = []; ballStickMap   = []; 
diffGif  = []; diffMap  = []; diffFigNo = 4;
montageGif  = []; montageMap  = []; montageFigNo = 5;
if ~isempty(pRes.priorMedFiltFigName)
    priorMedFiltFig = figure('visible', pRes.showFigs);
end
if ~isempty(pRes.ballStickGifName)
    ballStickFig = figure('Visible',pRes.showFigs);
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
rotAngleFromInd      = round(pRes.rotAngleFromInd, pRes.angleSigFig);
pRes.rotAngleFromInd = rotAngleFromInd;
assert(length(pRes.rotAngleFromInd) <= 255); % make sure that we can represent rotation with uint8

nBlocks = pRes.nBlockSpan^2;

blockZPriorFit = zeros(1,nBlocks);
blockZPriorData = zeros(1,nBlocks);
blockRPriorFit = zeros(1,nBlocks);
blockRPriorData = zeros(1,nBlocks);

% load stack (expect gif)
stackPath = fullfile(pRes.dataDir, subj, stackName);
moviePath = fullfile(pRes.dataDir, subj, movieName);
assert(exist(moviePath,'file') & exist(stackPath,'file'));

stackInf  = imfinfo(stackPath);
stackDepth = length(stackInf);

if ~isempty(pRes.loadedStack)
    stack = pRes.loadedStack;
else
    stack = zeros(stackInf(1).Height,stackInf(1).Width,stackDepth);
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
    movie = zeros(movieHeight,movieWidth,movieLength)
    for mm = 1:length(imfinfo(moviePath))
        movie(:,:,mm) = imread(moviePath,mm);
    end
end
rmfield(pRes,'loadedMovie');

% check stack and movie size matches file dimensions
assert(isequal([stackInf(1).Height, stackInf(1).Width, stackDepth], size(stack)));
assert(isequal([movieHeight, movieWidth, movieLength], size(movie)));

% crop dark edge off stack by removing rows of all zero entries and get new
% size
stack = cropStack(stack);
[stackHeight, stackWidth, stackDepth] = size(stack);

% create correlations directory
assert(moviePath(end-3)=='.');
movieDate = movieName(1:10);
movieDateDir = fullfile(pRes.dataDir, subj, movieDate);
corrDir = fullfile(movieDateDir, 'referenceLocalization');

if ~exist(movieDateDir, 'dir'), mkdir(movieDateDir); end
if ~exist(corrDir, 'dir'), mkdir(corrDir); end

if isempty(pRes.whichBlocks), pRes.whichBlocks = 1:nBlocks; end
if isempty(pRes.whichFrames), pRes.whichFrames = 1:movieLength; end
if isempty(pRes.whichSlices), pRes.whichSlices = 1:stackDepth; end


% divide movie into blocks based on the movie dimensions, desired number of
% blocks, percent overlap between blocks, and maximum amount of rotation
% (used to create a margin)
blockLocations = makeBlockLocations(movieHeight, movieWidth, ...
    pRes.nBlockSpan, pRes.blockOverlap, max(pRes.coarseRotAngles));

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
    frameCorrDir = fullfile(corrDir, sprintf('frame%03i', thisFrameNo));
    if ~exist(frameCorrDir, 'dir')
        mkdir(frameCorrDir);
    end
    
    
    [~, diffPlotMat]  = makeSubplots(diffFigNo,sqrt(nBlocks),sqrt(nBlocks),.1,.1,[ 0 0 1 1]);
    
    [~, montagePlotMat]  = makeSubplots(montageFigNo,sqrt(nBlocks),sqrt(nBlocks),.1,.1,[ 0 0 1 1]);
    
    
    if ff == 1
        xyzrPrior = getPrior(movieFrame, blockLocations, stack, pRes)
    else
        % might still want to have a filtering/interpolating step here in case a given frame goes haywire
        xyzrPrior = xyzrcPeak(1:5,:,thisFrameNo-1);
    end
    
    % ----------------- iterate through the blocks to keep values ------------------- %
    for bb = 1:length(pRes.whichBlocks)
        % select the relevant block location matrix and get its dimensions so
        % we can determine how many values will be in the correlation matrix
        thisBlockNo  = pRes.whichBlocks(bb);
        thisBlockLoc = blockLocations(:,:,thisBlockNo);
        bInf         = getBlockInf(thisBlockLoc);
        block        = movieFrame(bInf.indY,bInf.indX);
        
        % === Fill in correlations for neighborhood around peak ==== %
        % -----------------------------------------------------------%
        % specify a range of z values and rotation angles to save in the
        % persistent matrix
        thisNbrhdCtrZ = blockZPriorData(thisBlockNo);
        thisNbrhdCtrR = blockRPriorFit(thisBlockNo);
        
        zRangeToKeep = max(1,thisNbrhdCtrZ-ceil(pRes.nZToKeep/2)) : ...
            min(stackDepth,thisNbrhdCtrZ+ceil(pRes.nZToKeep/2));
        
        rotToKeep = fineRotAngles + thisNbrhdCtrR;
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
                    stackSlice, pRes.corrType));
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
        
        figure(11); clf
        plot(corrValsToSave); hold on
        scatter(peakInd, cPeak, 'mx');
        
        blockULInRefX = indX(peakInd) - bInf.width + 1;
        blkCtrInRefX  = double(blockULInRefX) + bInf.ctrX - 1;
        
        blockULInRefY = indY(peakInd) - bInf.height + 1;
        blkCtrInRefY  = double(blockULInRefY) + bInf.ctrY - 1;
        
        xyzrcPeak(:,thisBlockNo, thisFrameNo) = [blkCtrInRefX, blkCtrInRefY, ...
            double(indZ(peakInd)), double(indR(peakInd)),...
            double(cPeak)/double(intmax(pRes.corrType))]';
        
        
        
        % get the reference block from the peak location. this should
        % probably be its own function
        refBlockIndY = blockULInRefY : blockULInRefY + bInf.height - 1;
        refBlockIndX = blockULInRefX : blockULInRefX + bInf.width - 1;
        topPad = sum(refBlockIndY < 1);
        bottomPad = sum(refBlockIndY > stackHeight);
        leftPad = sum(refBlockIndX < 1);
        rightPad = sum(refBlockIndX > stackWidth);
        
        refBlockIndY(refBlockIndY < 1 | refBlockIndY > stackHeight) = [];
        refBlockIndX(refBlockIndX < 1 | refBlockIndX > stackWidth) = [];
        
        refBlock = stack(refBlockIndY, refBlockIndX, indZ(peakInd));
        refBlock = padarray(refBlock, [topPad leftPad], 0, 'pre');
        refBlock = padarray(refBlock, [bottomPad rightPad], 0, 'post');
        
        % get the subscript position for this block in the grid
        [blocki, blockj] = ind2sub([sqrt(nBlocks), sqrt(nBlocks)], thisBlockNo);
        
        % plot block-ref montage
        axes(montagePlotMat(blocki,blockj));
        imshowpair(blockRot, refBlock, 'montage');
        
        % plot block-ref diff
        axes(diffPlotMat(blocki,blockj));
        imagesc(blockRot - refBlock);
        colormap(bone);
        
        

        
        % % TODO: make two separate functions for pair plots and block
        % znbrhds
        % pair plots
        % %         createPairsPlot(pairsFigNo, peakInd, indX, indY, indZ, indR, ...
        % %             double(corrValsToSave)./double(intmax(pRes.corrType)),...
        % %             frameString, [.2 .9]);
        % %         [pairsGif, pairsMap]  = createGif(pairsFigNo, ff, movieLength, ...
        % %             pairsGif, pairsMap, fullfile(frameCorrDir, fprintf(pRes.pairsGifName,thisBlockNo)));
        %         % block Z neighbordhood report
        %     [~, blockZNbrhdPlotMat] = makeSubplots(blockZNbrhdFigNo, 4, 7, .1, .1, [ 0 0 1 1]);
        
        %         blockZNbrhdGif = []; blockZNbrhdMap = [];
        %         zNbrhdToPlot = indZ(peakInd)-3:indZ(peakInd)+3;
        %         for zz = 1:length(zNbrhdToPlot)
        %             thisZ = zNbrhdToPlot(zz);
        %             if thisZ < 1 | thisZ > stackDepth, continue; end
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
        %             blockZNbrhdGif,blockZNbrhdMap,fullfile(corrDir,blockZNbrhdGifName));
        
        % save block
        if ~isempty(pRes.blockSaveFormat)
            save(fullfile(frameCorrDir, sprintf(pRes.blockSaveFormat,thisBlockNo)), ...
                'corrValsToSave', 'indX', 'indY', 'indZ', 'indR', 'analysisDate', ...
                'stackPath','rotAngleFromInd','thisBlockLoc');
        end
        
        
        
    end
    
    % ball and stick image plot circles for the block centers with size
    %   proportional to rotation index and colored by their correlation value
    figure(ballStickFig); clf;
    thisFramePeaks = xyzrcPeak(:,:,thisFrameNo);
    scatter3(reshape(thisFramePeaks(1,:),[],1),reshape(thisFramePeaks(2,:),[],1),...
        reshape(thisFramePeaks(3,:),[],1), 400*(reshape(thisFramePeaks(4,:),[],1)),...%max(reshape(thisFramePeaks(4,:),[],1),1),...
        reshape(thisFramePeaks(5,:),[],1) , '.','linewidth',2);
    hold on
    % connect the circles with lines
    thisFramePeaksGridX = reshape(thisFramePeaks(1,:), sqrt(nBlocks), sqrt(nBlocks));
    thisFramePeaksGridY = reshape(thisFramePeaks(2,:), sqrt(nBlocks), sqrt(nBlocks));
    thisFramePeaksGridZ = reshape(thisFramePeaks(3,:), sqrt(nBlocks), sqrt(nBlocks));
    for bb = 1:sqrt(nBlocks)
        plot3(thisFramePeaksGridX(bb,:),thisFramePeaksGridY(bb,:),thisFramePeaksGridZ(bb,:),'-','color',[.6 .6 .5]);
        plot3(thisFramePeaksGridX(:,bb),thisFramePeaksGridY(:,bb),thisFramePeaksGridZ(:,bb),'-','color',[.6 .6 .5]);
    end
    xlabel('x');set(gca,'xdir','rev');xlim([0 1000]);ylim([0 1000])
    title(frameString);
    
    set(gca,'zdir','reverse','zlim',[1 stackDepth],'xlim',[0 stackWidth+bInf.width+1],...
        'ylim',[0 stackHeight+bInf.height+1],'clim',[.15 .85]);
    colormap(colormapRedBlue); colorbar;
    
    [ballStickFig, ballStickMap] = createGif(ballStickFigNo,ff,movieLength, ballStickFig,...
        ballStickMap, fullfile(corrDir,pRes.ballStickGifName));
    
    
    
    % turn figures into gifs
    %[rFitGif, rFitMap] = createGif(rFitFigNo,ff, movieLength,rFitGif,rFitMap,fullfile(corrDir,pRes.rFitGifName));
    %[zFitGif, zFitMap] = createGif(zFitFigNo,ff, movieLength,zFitGif,zFitMap,fullfile(corrDir, pRes.zFitGifName));
    [montageGif, montageMap] = createGif(montageFigNo,ff, movieLength,montageGif,montageMap,fullfile(corrDir, pRes.montageGifName));
    [diffGif, diffMap] = createGif(diffFigNo,ff, movieLength,diffGif,diffMap,fullfile(corrDir, pRes.diffGifName));
    
    
    % save summary
    if ~isempty(pRes.summarySaveName)
        save(fullfile(corrDir, pRes.summarySaveName), 'xyzrcPeak', 'blockLocations', ...
            'rotAngleFromInd','stackPath','moviePath','analysisDate');
    end
    
end

set(0, 'DefaultFigureVisible', 'on');

end


function [im, map]  = createGif(figNum, frameNo, movieLength, im, map, gifName)

% note on hardcopy: matlab says not to use it, but the alternative getframe
% requires you to have the figure open and fully visible to the user.
% the other alternative is print, but I can't figure how to collect
% its result in a variable
f.cdata = hardcopy(figure(figNum), '-Dzbuffer', '-r0');
if frameNo == 1
    [im, map] = rgb2ind(f.cdata,256,'nodither');
    im(1,1,1,movieLength) = 0;
else
    im(:,:,1,frameNo) = rgb2ind(f.cdata,map,'nodither');
end
if ~isempty(gifName)
    imwrite(im,map,gifName,'DelayTime',.2,'LoopCount',Inf);
end
end

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
nPeaks = length(peakInd);

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


function [xyzrPrior] = getPrior(movieFrame, blockLocations, whichBlocks, stack, pRes.whichSlices, pRes)
    
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

[stackHeight, stackWidth, stackDepth] = size(stack);
nBlocks = pRes.nBlockSpan^2;

% ------------ get initial rotation angle estimate for each block --------------- %
fprintf('getting initial prior on z for all blocks...\n');
bestZData = zeros(1,nBlocks);
bestZFit  = zeros(1,nBlocks);
bestRData = zeros(1,nBlocks);
bestRFit  = zeros(1,nBlocks);
bestXData = zeros(1,nBlocks);
bestYData = zeros(1,nBlocks);

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
    frameCorrVolNoRot = zeros(bInf.height+stackHeight-1,bInf.width+stackWidth-1,...
        stackDepth);
    
    for zz = 1:length(pRes.whichSlices),
        thisSliceNo = pRes.whichSlices(zz);
        stackSlice = stack(:,:,thisSliceNo);
        frameCorrVolNoRot(:,:,thisSliceNo) = computeBlockImageCorrs(block, stackSlice, 'double');
    end
    toc
    
    % find maximum correlation for each z and fit a polynomial to the
    % correlation values in a window around the max. Then choose z that
    % maximizes the fit curve
    maxCorrByZData = squeeze(max(max(frameCorrVolNoRot(:,:,:),[],1),[],2));
    [~, bestZData(thisBlockNo)] = max(maxCorrByZData);
    
    zIndToFit   = max(min(pRes.whichSlices), bestZData(thisBlockNo)-ceil(pRes.inferZWindow/2)) : ...
        min(max(pRes.whichSlices), bestZData(thisBlockNo)+ceil(pRes.inferZWindow/2));

    fitZ = polyfit(double(zIndToFit), double(maxCorrByZData(zIndToFit)'), pRes.zFitPower);
    
    maxCorrByZFit = polyval(fitZ, zIndToFit);
    bestZFit = zIndToFit(maxCorrByZFit == max(maxCorrByZFit));
    
    % save a report of the peak correlations and fits that will be used to get neighborhoods 
    if ~isempty(pRes.zFitFigName)
        [blocki, blockj] = ind2sub([pRes.nBlockSpan, pRes.nBlockSpan], thisBlockNo);
        
        ax1 = zFitPlotMat(blocki, blockj); hold(ax1, 'on');
        set(ax1,'box','off','ylim',[-.25 1], 'YTick',[max(maxCorrByZData)],...
            'xlim',[1 stackDepth], 'XTick', bestZData(thisBlockNo) , 'fontsize',6,...
            'xaxislocation','origin');                          
        
        % plot fit to data 
        plot(ax1,[bestZFit bestZFit],[0 1], '-.', ...
            zIndToFit, maxCorrByZFit, '-', 'color', [.9 .5 .7],'linewidth',1);
        % plot data points with mark best z
        plot(ax1,repmat(bestZData(thisBlockNo),1,2)],[0 1], '--','color', [.3 .5 .7],'linewidth',1);
        plot(ax1, 1:stackDepth, maxCorrByZData, '.', 'color', [.2 .4 .6], 'markersize',1.75)%;,...,
            %'markeredgecolor','k','linewidth',.05,'markersize',3);
        
        if blockj == 1 & blocki == ceil(pRes.nBlockSpan/2), ylabel(ax1,'correlation');
        end
        if blockj == ceil(pRes.nBlockSpan/2) & blocki == pRes.nBlockSpan, xlabel(ax1,'stack slice #')
        end
    end
   
    clear frameCorrVolNoRot;
end

if ~isempty(pRes.zFitFigName), saveas(zFitFig, fullfile(frameCorrDir, pRes.zFitFigName)); end

% smooth z estimate with median filter
zFit      = reshape(blockZPriorFit,sqrt(nBlocks),sqrt(nBlocks));
zData     = reshape(blockZPriorData,sqrt(nBlocks),sqrt(nBlocks));

zFitFilt  = medfilt2(zFit,[4 4], 'symmetric');
zDataFilt = medfilt2(zData,[4 4], 'symmetric');

% store z prior
if pRes.zPriorUseFit
    xyzrPrior(3,:) = round(reshape(zFitFilt,[],1));
else
    xyzrPrior(3,:) = reshape(zDataFilt,[],1));
end

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
    frameCorrVolRot = zeros(bInf.height+stackHeight-1,bInf.width+stackWidth-1,...
        length(pRes.coarseRotAngles), 'double');
    
    % get correlations for a pre-determined set of angles
    for rr = 1:length(pRes.coarseRotAngles)
        rotAngle = pRes.coarseRotAngles(rr);
        blockRot = rotateAndSelectBlock(movieFrame, thisBlockLoc, rotAngle);
        frameCorrVolRot(:,:,rr) = computeBlockImageCorrs(blockRot, ...
            bestStackSliceNoRot, 'double');
    end
    toc
    
    % find peak correlation for rotations using same process as for z
    [~, frameCorrVolMaxIx] = max(frameCorrVolRot(:));
    [bestYData(thisBlockNo), bestXData(thisBlockNo), bestRIndData] =  ind2sub(size(frameCorrVolRot), frameCorrVolMaxIx);
    
    bestYData(thisBlockNo) = bestYData(thisBlockNo) - bInf.height + bInf.ctrY;
    bestXData(thisBlockNo) = bestXData(thisBlockNo) - bInf.width + bInf.ctrX;

    bestRData(thisBlockNo) = pRes.coarseRotAngles(bestRIndData);
    maxCorrByRData = squeeze(max(max(frameCorrVolRot,[],1),[],2));
    
    fitRot =  polyfit((pRes.coarseRotAngles).',maxCorrByRData, pRes.rotFitPower);
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
        set(ax1,'box','off','ylim',[lowestCorr highestCorr], 'YTick',[max(maxCorrByRData)],...
            'xlim',[pRes.coarseRotAngles(1) pRes.coarseRotAngles(end)], 'XTick',[bestRData(thisBlockNo)],...
            'fontsize', 6, 'xaxislocation','origin'); 
        
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

if ~isempty(pRes.rFitFigName), saveas(rFitFig, fullfile(frameCorrDir, pRes.rFitFigName)); end


% smooth z estimate with median filter
rFit      = reshape(bestRFit,sqrt(nBlocks),sqrt(nBlocks));
rData     = reshape(blockRData,sqrt(nBlocks),sqrt(nBlocks));
rFitFilt  = medfilt2(rFit,[4 4], 'symmetric');
rDataFilt = medfilt2(rData,[4 4], 'symmetric');


% store z prior
if pRes.rPriorUseFit
    xyzrPrior(4,:) = reshape(rFitFilt,[],1);
else
    xyzrPrior(4,:) = reshape(rDataFilt,[],1));
end


xData = reshape(bestXData,pRes.nBlockSpan,pRes.nBlockSpan);
yData = reshape(bestYData,pRes.nblockSpan,pRes.nBlockSpan);
xDataFilt = medfilt2(xData,[4 4],'symmetric');
yDataFilt = medfilt2(yData,[4 4],'symmetric');

if ~isempty(pRes.priorMedFiltFigName)
    set(0, 'currentfigure', priorMedFiltFig);
    
    imagesc([zFit zFitFilt],'parent',subplot(2,3,1)); hold on;
    plot(repmat(mean(get(gca,'xlim')),1,2),get(gca,'ylim'),'k')
    caxis([min(pRes.whichSlices) max(pRes.whichSlices)]);  colorbar;
    title('Z Fits')
    
    imagesc([zData zDataFilt], 'parent',subplot(2,3,2)); hold on; 
    plot(repmat(mean(get(gca,'xlim')),1,2),get(gca,'ylim'),'k')
    caxis([min(pRes.whichSlices) max(pRes.whichSlices)]); colorbar;
    title('Z Data')
    
    imagesc([rFit rFitFilt],'parent',subplot(2,3,3)); hold on;
    plot(repmat(mean(get(gca,'xlim')),1,2),get(gca,'ylim'),'k')
    caxis([min(rFit(:)) max(rFit(:))]);%caxis([min(pRes.coarseRotAngles) max(pRes.coarseRotAngles)]);
    colorbar;
    title('R Fits')
    
    imagesc([rData rDataFilt], 'parent',subplot(2,3,4)); hold on;
    plot(repmat(mean(get(gca,'xlim')),1,2),get(gca,'ylim'),'k')
    %caxis([min(pRes.coarseRotAngles) max(pRes.coarseRotAngles)]);
    caxis([min(rFit(:)) max(rFit(:))]);
    colorbar;
    title('R Data')
    
    imagesc([yData yDataFilt],'parent',subplot(2,3,5)); hold on;
    plot(repmat(mean(get(gca,'xlim')),1,2),get(gca,'ylim'),'k')
    caxis([min(yData(:)) max(yData(:))]);%caxis([min(pRes.coarseRotAngles) max(pRes.coarseRotAngles)]);
    colorbar;
    title('Y Data')

    imagesc([xData xDataFilt],'parent',subplot(2,3,6)); hold on;
    plot(repmat(mean(get(gca,'xlim')),1,2),get(gca,'ylim'),'k')
    caxis([min(xData(:)) max(xData(:))]);%caxis([min(pRes.coarseRotAngles) max(pRes.coarseRotAngles)]);
    colorbar;
    title('X Data')
    
    colormap(redbluecmap);
    saveas(priorMedFiltFig, fullfile(corrDir, pRes.priorMedFiltFigName));
    
end

end

