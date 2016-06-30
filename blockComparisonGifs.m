function [montageGif, diffGif, overlapGif] = blockComparisonGifs(subj, movieDate)

useJukebox = exist('/jukebox','dir');
if useJukebox, 
    dataDir = '/jukebox/tank/jlgauthi/Data/';
else
    dataDir = '/Volumes/tank/jlgauthi/Data/';
end

% load peak locations and other useful information
movieDateDir = fullfile(dataDir, subj, movieDate);

corrDir = fullfile(movieDateDir, 'referenceLocalization');
load(fullfile(corrDir,'summary.mat'),'xyzrcoPeak','blockLocations','stackPath','pRes');
if ~exist('pRes','var'),
    load(fullfile(corrDir,'xyzrSearchRange.mat'),'pRes');
end
if useJukebox,
    stackPath = strrep(stackPath,'Volumes','jukebox');
else
    stackPath = strrep(stackPath,'jukebox','Volumes');
end

[~, nBlocks, nFrames] = size(xyzrcoPeak);



% load movie
moviePath    = sprintf('%s__L01__AVERAGE.tif',movieDateDir);
movieInf     = imfinfo(moviePath);
movie        = zeros(movieInf(1).Height,movieInf(1).Width,length(movieInf));
movieLength = size(movie,3);
for mm = 1:movieLength
    movie(:,:,mm) = imread(moviePath,mm);
end
% load stack
stackInf = imfinfo(stackPath);
stack = zeros(stackInf(1).Height,stackInf(1).Width,length(stackInf));
    for ss = 1:length(imfinfo(stackPath))
        stack(:,:,ss) = imread(stackPath,ss);
    end
stack = cropStack(stack);
[stackDim.height, stackDim.width, stackDim.depth] = size(stack);

% set up a grid for the montage figure and init some figures
maxBWidth = max(max(sum(blockLocations,2),[],1));
maxBHeight = max(max(sum(blockLocations,1),[],2));
montageGrid = zeros(( maxBHeight + 3) * pRes.nBlockSpan, ...
    2*pRes.nBlockSpan * (3+ maxBWidth));
montageFig = figure;
diffFig = figure;
overlapFig = figure;
montageGif = [];
montageMap = [];
diffGif = [];
diffMap = [];
overlapGif = [];
overlapMap = [];
[~, diffPlotMat]  = makeSubplots(diffFig,sqrt(nBlocks),sqrt(nBlocks),0,0,[ 0 0 1 1]);
[~, overlapPlotMat]  = makeSubplots(overlapFig,sqrt(nBlocks),sqrt(nBlocks),0,0,[ 0 0 1 1]);

frameOneHists = cell(nBlocks);

        

% normalize image versus reference
for thisFrameNo = 1:nFrames
    for thisBlockNo = 1:nBlocks
        [blocki, blockj] = ind2sub([sqrt(nBlocks), sqrt(nBlocks)], thisBlockNo);
        % get block info and block image
        movieFrame = movie(:,:,thisFrameNo);
        bInf = getBlockInf(blockLocations(:,:,thisBlockNo));
        bestBlockRot = rotateAndSelectBlock(movieFrame, bInf, xyzrcoPeak(4,thisBlockNo,thisFrameNo));
        
        % get the corresponding reference block iamge based on peak location
        blkCtrInRefX = xyzrcoPeak(1,thisBlockNo,thisFrameNo);
        blkCtrInRefY = xyzrcoPeak(2,thisBlockNo,thisFrameNo);
        refBlockIndX = blkCtrInRefX - (bInf.width-1)/2 : blkCtrInRefX + (bInf.width-1)/2;
        refBlockIndY = blkCtrInRefY - (bInf.height-1)/2 : blkCtrInRefY + (bInf.height-1)/2;
        % pad the reference block if it is smaller because the match partially overlaps with the edge of the reference image
        topPad = sum(refBlockIndY < 1);
        bottomPad = sum(refBlockIndY > stackDim.height);
        leftPad = sum(refBlockIndX < 1);
        rightPad = sum(refBlockIndX > stackDim.width);
        % zero out indices that don't land in the reference
        refBlockIndY(refBlockIndY < 1 | refBlockIndY > stackDim.height) = [];
        refBlockIndX(refBlockIndX < 1 | refBlockIndX > stackDim.width) = [];
        refBlock = stack(refBlockIndY, refBlockIndX, xyzrcoPeak(3,thisBlockNo,thisFrameNo));
        refBlock = padarray(refBlock, [topPad leftPad], 0, 'pre');
        refBlock = padarray(refBlock, [bottomPad rightPad], 0, 'post');
        
        % normalize blocks
        bestBlockRot = normalizeToZeroOne(bestBlockRot);
        refBlock     = normalizeToZeroOne(refBlock);
        
        % clip top and bottom 2 percent of image
        bestBlockRot(bestBlockRot>quantile(bestBlockRot(:),.98)) = quantile(bestBlockRot(:),.98);
        bestBlockRot(bestBlockRot<quantile(bestBlockRot(:),.02)) = quantile(bestBlockRot(:),.02);
        refBlock(refBlock>quantile(refBlock(:),.98)) = quantile(refBlock(:),.98);
        refBlock(refBlock<quantile(refBlock(:),.02)) = quantile(refBlock(:),.02);
        
        if thisFrameNo==1,
            frameOneHists{thisBlockNo} = bestBlockRot;
        else
            bestBlockRot    = imhistmatch(bestBlockRot,frameOneHists{thisBlockNo});
        end
        refBlock        = imhistmatch(refBlock,frameOneHists{thisBlockNo});
        %bestBlockRotFit = polyval(polyfit(bestBlockRot(:),refBlock(:),1),bestBlockRot);

        %         h=figure; imagesc([bestBlockRot refBlock bestBlockRotFit],'parent',h);
        %         close(h)
        % plot block-ref montage
        montageGrid(1+(maxBHeight+3)*(blocki-1):...
            (maxBHeight+3)*(blocki-1)+size(bestBlockRot,1),...
            1+(2*maxBWidth+3)*(blockj-1):...
            (2*maxBWidth+3)*(blockj-1)+2*size(bestBlockRot,2)) = ...
            [bestBlockRot refBlock];
        
        % montage image
        set(0,'currentfigure',montageFig);
        montageFigAx = cla;
        imagesc(montageGrid,'parent', montageFigAx);
        colormap(montageFigAx, 'bone');
        set(montageFig, 'Position', [1, 1, 1577, 954]);
        set(montageFigAx,'position',[.01 .01 .98 .98],'XTick',[],'YTick',[])
        axis(montageFigAx, 'image');

        
        
        % plot block-ref diff
        imagesc(bestBlockRot - refBlock,'parent',diffPlotMat(blocki,blockj));
        %imagesc(normalizeToZeroOne(bestBlockRot) - imhistmatch(normalizeToZeroOne(refBlock), normalizeToZeroOne(bestBlockRot)), 'parent',diffPlotMat(blocki,blockj));
        set(diffPlotMat(blocki,blockj), 'XTick',[],'YTick',[]);
        colormap(diffPlotMat(blocki,blockj), bone);
        axis(diffPlotMat(blocki,blockj),'image')
        % plot block-ref overlap
        imshowpair(bestBlockRot, refBlock,'falsecolor','parent',overlapPlotMat(blocki,blockj));
        axis(overlapPlotMat(blocki,blockj),'image')
        
        
        [diffGif, diffMap] = createGif(diffFig,thisFrameNo, movieLength,diffGif,...
            diffMap,fullfile(corrDir, 'blockDiffs.gif'));
        [montageGif, montageMap] = createGif(montageFig,thisFrameNo, movieLength,montageGif,...
            montageMap,fullfile(corrDir, 'blockMontage.gif'));
        [overlapGif, overlapMap] = createGif(overlapFig,thisFrameNo, movieLength,overlapGif,...
            overlapMap,fullfile(corrDir, 'blockOverlap.gif'));
    end
end
