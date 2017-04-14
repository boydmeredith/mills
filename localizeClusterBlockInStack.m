function [xyzrcoClusterPeaks, params] = ...
    localizeClusterBlockInStack(subject, theDate, location, varargin)
% identify where recorded ROIs are located in the reference stack
% 
% algorithm:
%
%   - use previous localization of the 1000-frame averages to approximate location
%   - perform localization of baseline average image, tightly centered on the approximate location
%   - compute transformation that maps every pixel in movie space to a pixel in reference space
%   - apply transformation to every ROI
%   - return transformed ROIs
%
%
%
% OUTPUTS:
%
%   roisTransformed - C x B struct array with fields
%                       w - Y x X x R matrix, weights of each ROI
%                       baseline - Y x X matrix, baseline image
%                       x - Y x X matrix, x-coordinate of each ROI pixel in the reference space
%                       y - Y x X matrix, y-coordinate of each ROI pixel in the reference space
%                       z - Y x X matrix, z-coordinate of each ROI pixel in the reference space
%                               this is the uniform within each block
%                       zFit - Y x X matrix, continuously-varying esimate of z-coordinate of each ROI pixel
%
%   xyzrcoClusterPeaks - 6 x C x B matrix, each entry the localization for a given cluster and block
%
%   roisRigid - same as roisTransformed, but with a limited rigid transformation (translate and rotate)
%
%   params - copy of parameters
%
%   extras - struct of extra things
%
%
%
% INPUTS: 
%
%   subject - e.g. 'J117'
%   theDate - e.g. '2015-12-06'
%   location - e.g. 'L01'
% 
%



p = inputParser;
addParamValue(p,'whichClusters',[]);
addParamValue(p,'whichBlocks',[]);
addParamValue(p,'xMargin',5);
addParamValue(p,'yMargin',5);
addParamValue(p,'zFitStyle','poly22');
addParamValue(p,'rFitStyle','poly11');
addParamValue(p,'zFitRobust','Bisquare');
addParamValue(p,'rFitRobust','Bisquare');

addParamValue(p,'loadedStackIs',[]);
addParamValue(p,'loadedStack',[]);
addParamValue(p,'refLocSumm',[]);

addParamValue(p,'zNbrhdRange',[-2:2]);
addParamValue(p,'rNbrhdRange',[-1:.25:1]);

% normalize stack brightness
addParamValue(p,'normalizeStack',true);

% normalize brightness of block image
addParamValue(p,'normalizeBlock',true);



addParamValue(p,'plot',[]);

%addParamValue(p,'xyzrcoClustPeaks',[]);
parse(p,varargin{:});

params = p.Results;



whichClusters = p.Results.whichClusters;
whichBlocks = p.Results.whichBlocks;
nbrhdInf.xMargin = p.Results.xMargin;
nbrhdInf.yMargin = p.Results.yMargin;
refLocSumm = p.Results.refLocSumm;
%stack = p.Results.stack;
nbrhdInf.zOffToKeep = p.Results.zNbrhdRange;
nbrhdInf.rOffToKeep = p.Results.rNbrhdRange;

meanCenter = @(x) x - mean(x);

% get some useful file names
nS = getNameStruct(subject, theDate, location);


% skip if clustering not done yet
if ~exist(nS.clusterFileName,'file')
    roisTransformed = [];
    roisRigid = [];
    xyzrcoClusterPeaks = [];
    return
end

% load the cluster file and figure out how many clusters there are
clusterFile = load(nS.clusterFileName,'clusterSpec');
nClusts = size(clusterFile.clusterSpec,2);

% load the first cluster info file and figure out how many blocks there are
for cc=1:nClusts
    if exist(nS.clusterInfoFileNameFcn(cc),'file')
        clustInfo = load(nS.clusterInfoFileNameFcn(cc));
        nBlocks = size(clustInfo.clusterBlockLocations,2);
        break
    end
end

% skip if no cluster info files were found
if ~exist('nBlocks','var')
    roisTransformed = [];
    roisRigid = [];
    xyzrcoClusterPeaks = [];
    return
end

% default to all clusters
if isempty(whichClusters)
    whichClusters = 1:nClusts;
else
    whichClusters = whichClusters(:)';
end

% default to all blocks
if isempty(whichBlocks)
    whichBlocks = 1:nBlocks;
else
    whichBlocks = whichBlocks(:)';
end

% load:
%   - localizations of each 1000-frame average
%   - path to the reference stack
%   - x-y locations of blocks of the 1000-frame averages
if isempty(refLocSumm)
    refLocSumm = load(nS.referenceLocalizationFileName,'xyzrcoPeak','stackPath','blockLocations');
end


% identify center of each block in the 1000-frame average space
blockCentersIn1000FrameAverage = nan(nBlocks,2);
for bb=1:nBlocks
    ctr = getBlockInf(refLocSumm.blockLocations(:,:,bb));
    blockCentersIn1000FrameAverage(bb,:) = [ctr.ctrX ctr.ctrY];
end



% load stack and get its dimensions
stackFileName = fullfile(PATH_DATA, refLocSumm.stackPath);
if ~isempty(p.Results.loadedStack)
    stack = p.Results.loadedStack;
else
    if isempty(p.Results.loadedStackIs)
        stackIs = imageSeries(stackFileName);
    else
        stackIs = p.Results.loadedStackIs;
    end
    stack = squeeze(permute(stackIs.images,[1 2 4 3]));
end


% initialize output variables
roisTransformed(nClusts, nBlocks) = struct('w',[],'x',[],'y',[],'z',[]);
roisRigid(nClusts, nBlocks) = struct('w',[],'x',[],'y',[],'z',[]);
xyzrcoClusterPeaks = nan(6, nClusts, nBlocks);


% store center points of all blocks
blockCentersX = nan(length(whichClusters),length(whichBlocks));
blockCentersY = nan(length(whichClusters),length(whichBlocks));


if p.Results.normalizeStack
    stackNorm = zeros(size(stack));
    for zz=1:size(stack,3)
        stackNorm(:,:,zz) = normalizeImageBrightness(double(stack(:,:,zz)));
    end
else
    stackNorm = stack;
end


for cc = whichClusters

    fprintf('\nworking on cluster %d, block...',cc);
    
    % list frames associated with this cluster
    thisClustFrames = find(matrixDownsampleSum(clusterFile.clusterSpec(:,cc),1000,1) > 0);
    
    if isempty(thisClustFrames), continue, end
    
    % get localization of 1000-frame averages corresponding to these frames
    x = squeeze(refLocSumm.xyzrcoPeak(1,:,thisClustFrames));
    y = squeeze(refLocSumm.xyzrcoPeak(2,:,thisClustFrames));
    z = squeeze(refLocSumm.xyzrcoPeak(3,:,thisClustFrames));
    r = squeeze(refLocSumm.xyzrcoPeak(4,:,thisClustFrames));
    
    % identify x-y offset between 1000-frame average space and the 512x512 space
    %   (load safe zone boundaries and use their initial values as offsets
    %   for the current block location centers)
    mc = load(nS.motCorrFileName);
    safeZoneX = computeMotionCorrectionBounds(mc.xShifts,[1 512],[1 512]);
    safeZoneY = computeMotionCorrectionBounds(mc.yShifts,[1 512],[1 512]);
    
    % add the offsets to the center points
    blockCentersIn512x512(:,1) = blockCentersIn1000FrameAverage(:,1) + safeZoneX(1);
    blockCentersIn512x512(:,2) = blockCentersIn1000FrameAverage(:,2) + safeZoneY(1);
    
    % fit a curve to the localization for these frames
    %   (this fit will be evaluated below to determine the center of a search range for final localization)
    xForFit = repmat(blockCentersIn512x512(:,1),[1 length(thisClustFrames)]);
    yForFit = repmat(blockCentersIn512x512(:,2),[1 length(thisClustFrames)]);
    fitForLeashCenter_x = fit([xForFit(:) yForFit(:)], x(:),'poly22');
    fitForLeashCenter_y = fit([xForFit(:) yForFit(:)], y(:),'poly22');
    fitForLeashCenter_z = fit([xForFit(:) yForFit(:)], z(:),'poly22');
    fitForLeashCenter_r = fit([xForFit(:) yForFit(:)], r(:),'poly22');
    
    
    
    for bb = whichBlocks
        if ~mod(bb,6)
            fprintf('  %d...',bb);
        end

        
        
        % load:
        %   x-y position of this block (in movie space)
        %   average baseline image 
        %   which 1000-frame averages are in this cluster
        [cframes, blockCenterX, blockCenterY, baselineImage]= ...
            getClusterBlockLocalization(subject,theDate,location,cc,bb);
        blockCentersX(cc,bb) = blockCenterX;
        blockCentersY(cc,bb) = blockCenterY;
        
        % skip if no baseline image (indicates baseline not computed)
        if isempty(baselineImage), continue, end
        
        % verify consitency
        assert(isequal(thisClustFrames,cframes));
        
        
        cImSz = size(baselineImage);
        
        % load the rois for this cluster and block
        try
            cellfile  = load(nS.cellFileNameFcn(cc,bb),'rois');
        catch
            warning(['Warning: Could not find'...
                'cluster file for cluster %03d, block %03d'],cc,bb);
            continue
        end
        
        % skip if no ROIs found
        if isempty(cellfile.rois), continue, end
        
        assert(isequal(cImSz, [size(cellfile.rois,1), size(cellfile.rois,2)]));
        
        
        
        % PERFORM FINAL LOCALIZATION OF THIS BLOCK IN REFERENCE SPACE
        
        
        % set up the neighborhood info struct to constrain possible localizations
        % evaluate the fit at the actual block center point
        nbrhdInf.xCtr = fitForLeashCenter_x(blockCenterX,blockCenterY);
        nbrhdInf.yCtr = fitForLeashCenter_y(blockCenterX,blockCenterY);
        nbrhdInf.zCtr = fitForLeashCenter_z(blockCenterX,blockCenterY);
        nbrhdInf.rCtr = fitForLeashCenter_r(blockCenterX,blockCenterY);
        
        % set up a binary block location matrix that indicates the subset
        % of the cluster block image that can stand the maximal rotation
        biggestAngle = max(abs(nbrhdInf.rCtr + nbrhdInf.rOffToKeep));
        rotatableSection = makeBlockLocations(cImSz(1),cImSz(2),[1 1], 0, biggestAngle);
        
        % normalize block image, if desired
        if p.Results.normalizeBlock
            cImNorm = normalizeImageBrightness(baselineImage);
        else
            cImNorm = baselineImage;
        end
        
        % localize this cluster block to the stack
        xyzrcoClusterPeaks(:,cc,bb) = localizeBlockInStackNbrhd(rotatableSection,...
            cImNorm, stackNorm, nbrhdInf);        
    end
end


% SAVE AUTOMATIC LOCALIZATION


% generate save name
filenameToSave = [nS.prefix '_referenceLocalizationBaseline.mat'];

% archive existing file, if there is one
archiveFile(filenameToSave,false)

% save
xyzrcoClusterPeaks_auto = xyzrcoClusterPeaks;
save(filenameToSave,'xyzrcoClusterPeaks_auto','params')


