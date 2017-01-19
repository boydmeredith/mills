function [roisTransformed, xyzrcoClusterPeaks, roisRigid, params, extras] = registerRois(subject, theDate, location, varargin)
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

%addParamValue(p,'stack',[]);
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

extras = struct;

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
stackIs = imageSeries(stackFileName);
stack = squeeze(permute(stackIs.images,[1 2 4 3]));


% initialize output variables
roisTransformed(nClusts, nBlocks) = struct('w',[],'x',[],'y',[],'z',[]);
roisRigid(nClusts, nBlocks) = struct('w',[],'x',[],'y',[],'z',[]);
xyzrcoClusterPeaks = nan(6, nClusts, nBlocks);


% store center points of all blocks
blockCentersX = nan(length(whichClusters),length(whichBlocks));
blockCentersY = nan(length(whichClusters),length(whichBlocks));


% store leash centers (for debugging)
leashCenters = nan(length(whichClusters),nBlocks,4);



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
        
        
        % note the leash center point (for debugging)
        leashCenters(cc,bb,:) = [nbrhdInf.xCtr nbrhdInf.yCtr nbrhdInf.zCtr nbrhdInf.rCtr];
        
        
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
        
        % extract values at peak correlation
        bestAngle = xyzrcoClusterPeaks(4,cc,bb);
        bestXCtr = xyzrcoClusterPeaks(1,cc,bb);
        bestYCtr = xyzrcoClusterPeaks(2,cc,bb);
        bestZ = xyzrcoClusterPeaks(3,cc,bb);
        
        
        
        
        
        % PERFORM RIGID TRANSFORMATION OF EACH ROI AND BASELINE IMAGE 
        
        
        % figure out how big cIm is going to be when rotated appropriately
        [rotH, rotW] = dimAfterRotation(cImSz(1), cImSz(2), bestAngle);
        padSz = ceil([rotH-cImSz(1)  rotW-cImSz(2)]./2);
        
        % pad roi weights with nans so they can rotate
        roiWPadded = padarray(cellfile.rois,padSz,nan,'both');
        roiWLoc = true(size(roiWPadded,1),size(roiWPadded,2));
        %nonNanLoc(~isnan(roiWPadded(:,:,1))) = true;
        
        % apply transformation
        roiWPaddedSz = size(roiWPadded);
        rotatedRoiW = nan(roiWPaddedSz);
        for ww = 1:roiWPaddedSz(3)
            rotatedRoiW(:,:,ww) = rotateAndSelectBlock(roiWPadded(:,:,ww),roiWLoc,bestAngle);
        end
        
        % transform baseline image
        cImPadded = padarray(baselineImage,padSz,nan,'both');
        blockBaselineTransformed = rotateAndSelectBlock(cImPadded,roiWLoc,bestAngle);
        
        % store in output variable
        roisRigid(cc,bb).w = rotatedRoiW;
        roisRigid(cc,bb).x = repmat(bestXCtr + meanCenter(1:roiWPaddedSz(2)),...
            roiWPaddedSz(1),1);
        roisRigid(cc,bb).y = repmat(bestYCtr + meanCenter(1:roiWPaddedSz(1))'...
            , 1,roiWPaddedSz(2));
        roisRigid(cc,bb).z =  repmat(bestZ,roiWPaddedSz(1),roiWPaddedSz(2));
        roisRigid(cc,bb).baseline = blockBaselineTransformed;
        
        
        
        
        
        % plot leash center vs final localization, if desired
        if ~isempty(p.Results.plot)
            figure(p.Results.plot);clf
            plot(leashCenters(cc,:,1),leashCenters(cc,:,2),'.')
            hold on;plot(squeeze(xyzrcoClusterPeaks(1,cc,:)),squeeze(xyzrcoClusterPeaks(2,cc,:)),'ro')
            axis ij
            drawnow
        end
        
    end
    
    
    
    
    
    
    
    
    
    % USE FINAL LOCALIZATION TO COMPUTE TRANSFORMATION FROM MOVIE SPACE TO REFERENCE SPACE
    
    
    % get summary of xyz points for this cluster
    blockInRefX = reshape(xyzrcoClusterPeaks(1,cc,:),[],1);
    blockInRefY = reshape(xyzrcoClusterPeaks(2,cc,:),[],1);
    blockInRefZ = reshape(xyzrcoClusterPeaks(3,cc,:),[],1);
    
    % choose only points that are not oddballs
    fitPoints = ~reshape(xyzrcoClusterPeaks(6, cc, :),[],1);
    
    % if insufficiently many non-oddballs, just use them all
    if sum(fitPoints) < nBlocks/2
       fitPoints = true(size(fitPoints));
    end
    
    % get transformation
    tform = cp2tform([ blockCentersX(cc,fitPoints)' blockCentersY(cc,fitPoints)'], [ blockInRefX(fitPoints) blockInRefY(fitPoints)], 'lwm');
    
    
    
    
    
    % APPLY TRANSFORMATION TO EACH BLOCK
    
    
    for bb = whichBlocks
        
        % get block range and average baseline image
        [~, ~, ~, baselineImage, blockRange]= getClusterBlockLocalization(subject,theDate,location,cc,bb);
        
        % transform baseline image and each ROI (all at once)
        cellfile  = load(nS.cellFileNameFcn(cc,bb),'rois');
        [dataTransformed, xData, yData] = imtransform(cat(3,baselineImage,cellfile.rois),tform,...
            'udata',blockRange.blockRangeX,'vdata',blockRange.blockRangeY,'fillvalues', nan);
        
        bsImTransformed = dataTransformed(:,:,1);
        roiShapesTransformed = dataTransformed(:,:,2:end);
        
        % note the x and y values of pixels in the transformed image
        xVals = linspace(round(xData(1)),round(xData(2)),size(bsImTransformed,2));
        yVals = linspace(round(yData(1)),round(yData(2)),size(bsImTransformed,1));
        
        % set to nan where there is no image data
        xInd = repmat(xVals,length(yVals),1);
        yInd = repmat(yVals',1,length(xVals));
        xInd(isnan(bsImTransformed)) = nan;
        yInd(isnan(bsImTransformed)) = nan;
        
        roisTransformed(cc,bb).w = roiShapesTransformed;
        roisTransformed(cc,bb).x = xInd;
        roisTransformed(cc,bb).y = yInd;
        % use same Z value as above
        roisTransformed(cc,bb).z =  repmat(round(mean(roisRigid(cc,bb).z(:))),length(yVals),length(xVals));
        roisTransformed(cc,bb).baseline = bsImTransformed;
    end
    

    
    
    % fit Z (in reference frame space) so that Z values can be interpolated everywhere
    
    % if there are sufficiently many non-nan values...
    usePts = ~any(isnan([blockInRefX blockInRefY blockInRefZ]),2);
    if sum(usePts) > nBlocks/2
        
        % compute z plane fit
        fitZ = fit([blockInRefX(usePts) blockInRefY(usePts) ], blockInRefZ(usePts),...
            p.Results.zFitStyle,'Robust',p.Results.zFitRobust);
        
        % store result 
        for bb=whichBlocks
            roisRigid(cc,bb).zFit = fitZ(roisRigid(cc,bb).x,roisRigid(cc,bb).y);
            roisTransformed(cc,bb).zFit = fitZ(roisTransformed(cc,bb).x,roisTransformed(cc,bb).y);
        end
    else
        % otherwise just use the fit z values
        roisRigid(cc,bb).zFit = roisRigid(cc,bb).z;
        roisTransformed(cc,bb).zFit = roisTransformed(cc,bb).z;
    end
    
    disp(' ')
end


% populate extras
extras.leashCenters = leashCenters;
extras.blockCentersX = blockCentersX;
extras.blockCentersY = blockCentersY;
extras.tformMovieToReference = tform;
extras.fitZfromXY = fitZ;


if 0 % for debugging
    
    %
    figure;hist(leashCenters(cc,:,2)' - squeeze(xyzrcoClusterPeaks(2,cc,:)),50)
    
    %
    
    stackNorm = zeros(size(stack));
    for zz=1:size(stack,3)
        stackNorm(:,:,zz) = normalizeImageBrightness(double(stack(:,:,zz)));
    end
    %
    figure(101);clf
    plot(leashCenters(cc,:,1),leashCenters(cc,:,2),'.')
    hold on;plot(squeeze(xyzrcoClusterPeaks(1,cc,:)),squeeze(xyzrcoClusterPeaks(2,cc,:)),'ro')
    axis image
    
    %
    
    
    %
    
    
    
    figure;
    plot(cxList(cc,:),cyList(cc,:),'.')

    
    xErrors = leashCenters(cc,:,1)' - blockInRefX;
    yErrors = leashCenters(cc,:,2)' - blockInRefY;
    zErrors = leashCenters(cc,:,3)' - blockInRefZ;
    
    figure;
    plot(abs(zErrors), abs(xErrors),'o')
    
    text(abs(zErrors), abs(xErrors),strsplit(num2str(1:36))')
    ylabel('x error')
    xlabel('z error')
end
