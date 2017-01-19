function [M, C] = loadCells(mouse, theDate, clusterNo, location)
% function [M, C] = loadCells(mouse, theDate,clusterNo, location)

% input:
% mouse, 
%
% output:
% M     is a matrix with size 512 x 512 x nCells. It contains 512 x 512 images
%       of each cell
%
% C     is a matrix with size 512 x 512 x nCells. It contains 512 x 512 images 
%       of each cell where only the centroid is labeled with the cellId


% load the name struct
nS = getNameStruct(mouse,theDate,location);

% load the cellMatchMatrix
cellMatchFile = fullfile(PATH_DATA,mouse,'cellMatchMatrix.mat');
cellMatch = load(cellMatchFile);

% load the clusterInfoFile
clustInfo = load(nS.clusterInfoFileNameFcn(clusterNo));

% find the number that corresponds to this date
dateNo = find(strcmp(cellMatch.dateList,theDate));

sessionInd = (cellMatch.sessionNumbers(:,1) == dateNo & ...
    cellMatch.sessionNumbers(:,2) == clusterNo);


% cellMatchMatrix is nCells x nSess x 2
% the last dimension specifies the block number and the cell id within that
% block
blockNumbers = unique(cellMatch.cellMatchMatrix(:,sessionInd,1));
blockNumbers(blockNumbers==0) = [];

% load reference localization summary file
refLocSum = load(nS.referenceLocalizationFileName);

% for each block and for each cell
M = nan(512,512, 1000);
C = nan(512,512, 1000);

for bb = blockNumbers'
    
    % determine which entries in cellMatchMatrix correspond to this session and block
    blockInd = cellMatch.cellMatchMatrix(:,sessionInd,1) == bb;
    
    % determine within block cell ids based on entry in cellMatchMatrix(:,:,2)
    cellsInBlock = cellMatch.cellMatchMatrix(blockInd,sessionInd,2);
    
    % determine global cell ids based on index into the cell match matrix
    cellId = find(blockInd);
    
    % load the images of the rois for this block    
    cellFile = load(nS.cellFileNameFcn(clusterNo,bb));
    thisRois = cellFile.rois;
    
    % determine the x and y coordinates of this block
        % load the localized coordinates of this block in the volume
        thisBlockRangeY = clustInfo.clusterBlockLocations{bb}.blockRangeY;
        thisBlockRangeX = clustInfo.clusterBlockLocations{bb}.blockRangeX;
        yInd = thisBlockRangeY(1):thisBlockRangeY(2);
        xInd = thisBlockRangeX(1):thisBlockRangeX(2);

    
    % loop over the cells and put them into the roi matrix
    for ii = 1:length(cellsInBlock)
        thisGlobalId  = cellId(ii);
        thisIdInBlock = cellsInBlock(ii); 
        
        M(yInd, xInd, thisGlobalId) = thisRois(:,:,thisIdInBlock);
        binaryRoi = M(:,:, cellId)>0;
        rp = regionprops(binaryRoi);
        C(round(rp.Centroid(2)),round(rp.Centroid(1)),thisGlobalId) = thisGlobalId;
        
    end
        
end


% TO DO:
%
% images of the blocks live in the baseline file
% bsFile.extraStats.F_mean is the image of the rois







