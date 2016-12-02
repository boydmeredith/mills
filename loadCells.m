function M = loadCells(mouse, theDate,theSessionNo, location)
% making this function so that we can load all rois from a day into image
% frame coordinates. we will use this to do image frame registration with
% turbo reg. it would be easy to add functionality to load into reference
% space instead of image frame space.


% load the name struct
nS = getNameStruct(mouse,theDate,location);

% load the cellMatchMatrix
cellMatchFile = fullfile(PATH_DATA,mouse,'cellMatchMatrix.mat');
cellMatch = load(cellMatchFile);

% load the clusterInfoFile
clustInfo = load(nS.clusterFileName);

dateIndInDateList = find(strcmp(cellMatch.dateList,theDate));

sessionInd = (cellMatch.sessionNumbers(:,1) == dateIndInDateList & ...
    cellMatch.sessionNumbers(:,2) == theSessionNo);


% cellMatchMatrix is nCells x nSess x 2
% the last dimension specifies the block number and the cell id within that
% block
blockNumbers = unique(cellMatch.cellMatchMatrix(:,sessionInd,1));
blockNumbers(blockNumbers==0) = [];


% for each block and for each cell

for bb = blockNumbers'
    blockInd = cellMatch.cellMatchMatrix(:,sessionInd,1) == bb;
    cellsInBlock = cellMatch.cellMatchMatrix(blockInd,sessionInd,2);
    
end
% find the position of the block and put the cell into the 512x512 coords
% of the full image









