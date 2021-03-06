function [isSt] = blockComparisonGifs(subj, movieDate, location, whichCompsToSave)
% function [isSt] = blockComparisonGifs(subj, movieDate, location, whichCompsToSave)
%
% Create imageSeries showing comparisons of the blocks in a movie and the
% matches found in a reference. Image histograms are modified to get rid of
% uninformative flickering from frame to frame.
%
% Input: 
% - subj, movieDate, location
% - whichCompsToSave:   any combination of the letters 'mod'?they indicate
%                       which image series to save as gifs.
%
% Output: 
% Three comparisons are packaged into the struct isSt, they are:
% - overlap:    shows the reference in the Red channel and the blocks in the
%              (greenish)-Blue channel
% - montage:    shows the images side by side
% - diff:       shows the block - the found match 
%
%
if isempty(location), location = 'L01'; end
if nargin < 4, whichCompsToSave = ''; end
% load peak locations and other useful information
movieDateDir = fullfile(jlgDataDir, subj, movieDate);

corrDir = referenceLocalizationDir(subj, movieDate, location);
s       = load(fullfile(corrDir,'summary.mat'),'xyzrcoPeak','blockLocations',...
                                                    'stackPath','params');

fullStackPath = fullfile(jlgDataDir, s.stackPath);

[~, nBlocks, nFrames] = size(s.xyzrcoPeak);

% load movie
moviePath    = sprintf('%s__L01__AVERAGE.tif',movieDateDir);
movieInf     = imfinfo(moviePath);
movie        = zeros(movieInf(1).Height,movieInf(1).Width,length(movieInf));
movieLength = size(movie,3);
fprintf('loading movie for subj: %s date: %s location: %s\n', subj, movieDate,location);
for mm = 1:movieLength
    movie(:,:,mm) = imread(moviePath,mm);
end
% load stack
stackInf = imfinfo(fullStackPath);
stack = zeros(stackInf(1).Height,stackInf(1).Width,length(stackInf));
fprintf('loading stack...\n');
for ss = 1:length(stackInf)
    stack(:,:,ss) = imread(fullStackPath,ss);
end
stack = cropStack(stack);
[stackDim.height, stackDim.width, stackDim.depth] = size(stack);

% set up a grid for the montage figure and init some figures
maxBWidth = max(max(sum(s.blockLocations,2),[],1));
maxBHeight = max(max(sum(s.blockLocations,1),[],2));
montageGrid = zeros(5+ ( maxBHeight + 3) * s.params.mByNBlocks(1), ...
    2*s.params.mByNBlocks(2) * (3+ maxBWidth), movieLength);
diffGrid = zeros(5+ ( maxBHeight + 3) * s.params.mByNBlocks(1), ...
    s.params.mByNBlocks(2) * (3+ maxBWidth), movieLength);
overlapGrid = zeros(5+ ( maxBHeight + 3) * s.params.mByNBlocks(1), ...
    s.params.mByNBlocks(2) * (3+ maxBWidth), 3, movieLength);


% matrix to store histograms for the first frame for normalization of
% subsequent frames to avoid uninformative flickering
frameOneHists = cell(nBlocks);


for thisFrameNo = 1:nFrames
    fprintf('frame: %03i/%03i...\t',thisFrameNo,nFrames)
    for thisBlockNo = 1:nBlocks
        [blocki, blockj] = ind2sub(s.params.mByNBlocks, thisBlockNo);

        % get block info and block image
        movieFrame = movie(:,:,thisFrameNo);
        bInf = getBlockInf(s.blockLocations(:,:,thisBlockNo));
        bestBlockRot = rotateAndSelectBlock(movieFrame, bInf, s.xyzrcoPeak(4,thisBlockNo,thisFrameNo));
        
        % get the corresponding reference block iamge based on peak location
        blkCtrInRefX = s.xyzrcoPeak(1,thisBlockNo,thisFrameNo);
        blkCtrInRefY = s.xyzrcoPeak(2,thisBlockNo,thisFrameNo);
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
        refBlock = stack(refBlockIndY, refBlockIndX, s.xyzrcoPeak(3,thisBlockNo,thisFrameNo));
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

        % plot block-ref montage
        
        % set up the indices for the plots
        xStartMontage = 1+(2*maxBWidth+3)*(blockj-1);
        xEndMontage = xStartMontage+2*size(bestBlockRot,2)-1;
        
        yStart = 1+(maxBHeight+3)*(blocki-1);
        xStart = 1+(maxBWidth+3)*(blockj-1);
        
        yEnd   = yStart+size(bestBlockRot,1)-1;
        xEnd   = xStart+size(bestBlockRot,2)-1;
        
        
        montageGrid(yStart : yEnd, xStartMontage:xEndMontage, ...
            thisFrameNo) = [bestBlockRot refBlock];
        
        % plot diff image
        %diffGrid(yStart:yEnd,xStart:xEnd,thisFrameNo) = (bestBlockRot - refBlock)./max(abs((bestBlockRot(:) - refBlock(:))));
        diffGrid(yStart:yEnd,xStart:xEnd,thisFrameNo) = (bestBlockRot - refBlock)./max(abs(refBlock(:)));
        
        % plot falsecolor image
        overlapGrid(yStart:yEnd,xStart:xEnd,:,thisFrameNo) = cat(3,refBlock,bestBlockRot, bestBlockRot);
        
    end
    timePoint = round(thisFrameNo/nFrames*size(diffGrid,2));
    montageGrid(end-3:end,1:timePoint,thisFrameNo) = 1;
    diffGrid(end-3:end,1:timePoint,thisFrameNo) = 1;
    overlapGrid(end-3:end,1:timePoint,:,thisFrameNo) = 1;
    assert(max(diffGrid(:))<=1 && min(diffGrid(:))>=-1);
end
diffGrid(end,end,:)=-1;

isSt.overlap = imageSeries(overlapGrid);
isSt.diff    = imageSeries(diffGrid);
isSt.montage = imageSeries(montageGrid);
if ismember('o',whichCompsToSave)
    isSt.overlap.saveImages(fullfile(corrDir,'blockOverlap'),'imageType','gif','class','uint8','normalizeColors',true,'overwrite',true);
end
if ismember('m',whichCompsToSave)
    isSt.montage.saveImages(fullfile(corrDir,'blockmontage'),'imageType','gif','class','uint8','normalizeColors',true,'overwrite',true);
end
if ismember('d',whichCompsToSave)
    isSt.diff.saveImages(fullfile(corrDir,'blockDiff'),'imageType','gif','class','uint8','normalizeColors',true,'overwrite',true);
end
