function [roisRigid, roisNonrigid, xyzrcoClusterPeaks, params] = ...
    registerRois(subject, theDate, location, varargin)



p = inputParser;
addParamValue(p,'whichClusters',[]);
addParamValue(p,'whichBlocks',[]);
addParamValue(p,'xMargin',5);
addParamValue(p,'yMargin',5);
addParamValue(p,'zFitStyle','poly11');
addParamValue(p,'rFitStyle','poly11');
addParamValue(p,'zFitRobust','Bisquare');
addParamValue(p,'rFitRobust','Bisquare');

%addParamValue(p,'stack',[]);
addParamValue(p,'refLocSumm',[]);

addParamValue(p,'zNbrhdRange',[-2:2]);
addParamValue(p,'rNbrhdRange',[-1:.25:1]);

% normalize stack brightness
addParamValue(p,'normalizeStack',false);

% normalize brightness of block image
addParamValue(p,'normalizeBlock',false);


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
    roisRigid = [];
    xyzrcoClusterPeaks = [];
    return
end

% load first cell finding file to get estimate for cell counts and block
% size
%cellfile = load(nS.cellFileNameFcn(1,1),'rois');

if isempty(whichClusters)
    whichClusters = 1:nClusts;
end
if isempty(whichBlocks)
    whichBlocks = 1:nBlocks;
end

% load registered peaks and the path to the stack used 
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
stackSz = size(stack);

% initialize output variables
xyzrcoClusterPeaks = nan(6, nClusts, nBlocks);
roisRigid(nClusts, nBlocks) = struct('w',[],'x',[],'y',[],'z',[]);
roisNonrigid(nClusts, nBlocks) = struct('w',[],'x',[],'y',[],'z',[]);


% store values for debuggin
cxList = nan(max(whichClusters),max(whichBlocks));
cyList = nan(max(whichClusters),max(whichBlocks));
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
    
    % find localization for these frames
    x=squeeze(refLocSumm.xyzrcoPeak(1,:,thisClustFrames));
    y=squeeze(refLocSumm.xyzrcoPeak(2,:,thisClustFrames));
    z=squeeze(refLocSumm.xyzrcoPeak(3,:,thisClustFrames));
    r=squeeze(refLocSumm.xyzrcoPeak(4,:,thisClustFrames));
    
    % identify offset between 1000-frame average space and the 512x512 space
    %   (load safe zone boundaries and use their initial values as offsets
    %   for the current block location centers)
    mc = load(nS.motCorrFileName);
    safeZoneX = computeMotionCorrectionBounds(mc.xShifts,[1 512],[1 512]);
    safeZoneY = computeMotionCorrectionBounds(mc.yShifts,[1 512],[1 512]);
    
    % add the offsets to the center points
    blockCentersIn512x512(:,1) = blockCentersIn1000FrameAverage(:,1) + safeZoneX(1);
    blockCentersIn512x512(:,2) = blockCentersIn1000FrameAverage(:,2) + safeZoneY(1);
    
    % fit a curve to the localization for these frames
    xForFit = repmat(blockCentersIn512x512(:,1),[1 length(thisClustFrames)]);
    yForFit = repmat(blockCentersIn512x512(:,2),[1 length(thisClustFrames)]);
    fitForLeashCenter_x = fit([xForFit(:) yForFit(:)], x(:),'poly22');
    fitForLeashCenter_y = fit([xForFit(:) yForFit(:)], y(:),'poly22');
    fitForLeashCenter_z = fit([xForFit(:) yForFit(:)], z(:),'poly22');
    fitForLeashCenter_r = fit([xForFit(:) yForFit(:)], r(:),'poly22');
    
    
    
    fitZ = fit([x(:) y(:)], z(:),p.Results.zFitStyle,'Robust',p.Results.zFitRobust);
    fitR = fit([x(:) y(:)], r(:),p.Results.rFitStyle,'Robust',p.Results.rFitRobust);
    
    
    
    for bb = whichBlocks
        if ~mod(bb,6)
            fprintf('  %d...',bb);
        end

        
        % get positioning and average image for block in the cluster
        [cframes, cx, cy, cIm]= getClusterBlockLocalization(subject,theDate,...
            location,cc,bb);
        cxList(cc,bb) = cx;
        cyList(cc,bb) = cy;
        
        if isempty(cIm), continue, end
        
        
        assert(isequal(thisClustFrames,cframes));
        
        cImSz = size(cIm);
        
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
        
        % set up the neighborhood info struct to constrain possible localizations
        nbrhdInf.xCtr = fitForLeashCenter_x(cx,cy);
        nbrhdInf.yCtr = fitForLeashCenter_y(cx,cy);
        nbrhdInf.zCtr = fitForLeashCenter_z(cx,cy);
        nbrhdInf.rCtr = fitForLeashCenter_r(cx,cy);
        
        leashCenters(cc,bb,:) = [nbrhdInf.xCtr nbrhdInf.yCtr nbrhdInf.zCtr nbrhdInf.rCtr];
        
        % delete!
        %nbrhdInf.xCtr = cx;
        %nbrhdInf.yCtr = cy;
        %nbrhdInf.zCtr = round(fitZ(cx,cy));
        %nbrhdInf.rCtr = fitR(cx,cy);
        
        
        
        % set up a binary block location matrix that indicates the subset
        % of the cluster block image that can stand the maximal rotation
        biggestAngle = max(abs(nbrhdInf.rCtr + nbrhdInf.rOffToKeep));
        rotatableSection = makeBlockLocations(cImSz(1),cImSz(2),[1 1], 0, biggestAngle);
        
        % normalize block image, if desired
        if p.Results.normalizeBlock
            cImNorm = normalizeImageBrightness(cIm);
        else
            cImNorm = cIm;
        end
        
        % localize this cluster block to the stack
        xyzrcoClusterPeaks(:,cc,bb) = localizeBlockInStackNbrhd(rotatableSection,...
            cImNorm, stackNorm, nbrhdInf);
        
        % extract values at peak correlation
        bestAngle = xyzrcoClusterPeaks(4,cc,bb);
        bestXCtr = xyzrcoClusterPeaks(1,cc,bb);
        bestYCtr = xyzrcoClusterPeaks(2,cc,bb);
        bestZ = xyzrcoClusterPeaks(3,cc,bb);
        
        % figure out how big cIm is going to be when rotated appropriately
        [rotH, rotW] = dimAfterRotation(cImSz(1), cImSz(2), bestAngle);
        padSz = ceil([rotH-cImSz(1)  rotW-cImSz(2)]./2);
        
        % pad roi weights with nans so they can rotate
        roiWPadded = padarray(cellfile.rois,padSz,nan,'both');
        roiWLoc = true(size(roiWPadded,1),size(roiWPadded,2));
        %nonNanLoc(~isnan(roiWPadded(:,:,1))) = true;
        
        
        % transform each ROI
        roiWPaddedSz = size(roiWPadded);
        rotatedRoiW = nan(roiWPaddedSz);
        for ww = 1:roiWPaddedSz(3)
            rotatedRoiW(:,:,ww) = rotateAndSelectBlock(roiWPadded(:,:,ww),...
                roiWLoc,bestAngle);
        end
        
        % transform baseline image
        cImPadded = padarray(cIm,padSz,nan,'both');
        blockBaselineTransformed = rotateAndSelectBlock(cImPadded,roiWLoc,bestAngle);
        
        % store in output variable
        roisRigid(cc,bb).w = rotatedRoiW;
        roisRigid(cc,bb).x = repmat(bestXCtr + meanCenter(1:roiWPaddedSz(2)),...
            roiWPaddedSz(1),1);
        roisRigid(cc,bb).y = repmat(bestYCtr + meanCenter(1:roiWPaddedSz(1))'...
            , 1,roiWPaddedSz(2));
        roisRigid(cc,bb).z =  repmat(bestZ,roiWPaddedSz(1),roiWPaddedSz(2));
        roisRigid(cc,bb).baseline = blockBaselineTransformed;
        
        % plot, if desired
        if ~isempty(p.Results.plot)
            figure(p.Results.plot);clf
            plot(leashCenters(cc,:,1),leashCenters(cc,:,2),'.')
            hold on;plot(squeeze(xyzrcoClusterPeaks(1,cc,:)),squeeze(xyzrcoClusterPeaks(2,cc,:)),'ro')
            axis ij
            drawnow
        end
        
    end
    
    % interpolate z value at each pixel by fitting Z
    xThisCluster = reshape(xyzrcoClusterPeaks(1,cc,:),[],1);
    yThisCluster = reshape(xyzrcoClusterPeaks(2,cc,:),[],1);
    zThisCluster = reshape(xyzrcoClusterPeaks(3,cc,:),[],1);
    fitZ = fit([xThisCluster yThisCluster ], zThisCluster,p.Results.zFitStyle,'Robust',p.Results.zFitRobust);
    
    % store result
    for bb=whichBlocks
        roisRigid(cc,bb).zFit = fitZ(roisRigid(cc,bb).x,roisRigid(cc,bb).y);
    end
    

    % interpolate z value to every point in the image space
    [xx,yy] = meshgrid(1:512);
    zFitInMovieSpace = fit([cxList' cyList'], zThisCluster,'poly22');
    zMovieSpace = zFitInMovieSpace(xx,yy);
    
    
    % compute overall transformation from image space to reference space
    goodOnes = ~reshape(xyzrcoClusterPeaks(6, cc, :),[],1);
    tform = cp2tform([ cxList(goodOnes)' cyList(goodOnes)'], [ xThisCluster(goodOnes) yThisCluster(goodOnes)], 'lwm');
    
    % for each block
    for bb = whichBlocks
        
        % get block range and average baseline image
        [~, ~, ~, cIm, blockRange]= getClusterBlockLocalization(subject,theDate,location,cc,bb);
        
        % transform baseline image
        [cImTransformed, xData, yData] = imtransform(cIm,tform,'udata',...
            blockRange.blockRangeX,'vdata',blockRange.blockRangeY,'fillvalues', nan);
        
        % transform each ROI
        cellfile  = load(nS.cellFileNameFcn(cc,bb),'rois');
        roisTransformed = imtransform(cellfile.rois,tform,'udata',...
            blockRange.blockRangeX,'vdata',blockRange.blockRangeY,'fillvalues', nan);
        
        
        % fit z values
        yPix = blockRange.blockRangeY(1):blockRange.blockRangeY(2);
        xPix = blockRange.blockRangeX(1):blockRange.blockRangeX(2);
        zFitValues = imtransform(zMovieSpace(yPix,xPix),tform,'udata',...
            blockRange.blockRangeX,'vdata',blockRange.blockRangeY,'fillvalues', nan);
        
        % note the x and y values of pixels in the transformed image
        xVals = linspace(round(xData(1)),round(xData(2)),size(cImTransformed,2));
        yVals = linspace(round(yData(1)),round(yData(2)),size(cImTransformed,1));
        
        roisNonrigid(cc,bb).w = roisTransformed;
        roisNonrigid(cc,bb).x = repmat(xVals,length(yVals),1);
        roisNonrigid(cc,bb).y = repmat(yVals',1,length(xVals));
        % use same Z value as above
        roisNonrigid(cc,bb).z =  repmat(round(mean(roisRigid(cc,bb).z(:))),length(yVals),length(xVals));
        roisNonrigid(cc,bb).zFit =  zFitValues;
        roisNonrigid(cc,bb).baseline = cImTransformed;
    end
    
end





%%
if 0 % for debugging
    
    %%
    figure;hist(leashCenters(cc,:,2)' - squeeze(xyzrcoClusterPeaks(2,cc,:)),50)
    
    %%
    
    stackNorm = zeros(size(stack));
    for zz=1:size(stack,3)
        stackNorm(:,:,zz) = normalizeImageBrightness(double(stack(:,:,zz)));
    end
    %%
    figure(101);clf
    plot(leashCenters(cc,:,1),leashCenters(cc,:,2),'.')
    hold on;plot(squeeze(xyzrcoClusterPeaks(1,cc,:)),squeeze(xyzrcoClusterPeaks(2,cc,:)),'ro')
    axis image
    
    %%
    
    
    %%
    
    
    
    figure;
    plot(cxList(cc,:),cyList(cc,:),'.')

    
    xErrors = leashCenters(cc,:,1)' - xThisCluster;
    yErrors = leashCenters(cc,:,2)' - yThisCluster;
    zErrors = leashCenters(cc,:,3)' - zThisCluster;
    
    figure;
    plot(abs(zErrors), abs(xErrors),'o')
    
    text(abs(zErrors), abs(xErrors),strsplit(num2str(1:36))')
    ylabel('x error')
    xlabel('z error')
end
