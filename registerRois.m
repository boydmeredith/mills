function [rois, xyzrcoClusterPeaks] = ...
    registerRois(subject, theDate, location, varargin)


p = inputParser;
addOptional(p,'whichClusters',[]);
addOptional(p,'whichBlocks',[]);
addOptional(p,'xMargin',15);
addOptional(p,'yMargin',15);
addOptional(p,'zFitStyle','poly11');
addOptional(p,'rFitStyle','poly11');
addOptional(p,'zFitRobust','Bisquare');
addOptional(p,'rFitRobust','Bisquare');

%addOptional(p,'stack',[]);
addOptional(p,'refLocSumm',[]);

addOptional(p,'zNbrhdRange',[-2:2]);
addOptional(p,'rNbrhdRange',[-1:.25:1]);

%addOptional(p,'xyzrcoClustPeaks',[]);
parse(p,varargin{:});

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

% load the cluster file and figure out how many clusters there are
clusterFile = load(nS.clusterFileName,'clusterSpec');
nClusts = size(clusterFile.clusterSpec,2);

% load the first cluster info file and figure out how many blocks there are
clustInfo = load(nS.clusterInfoFileNameFcn(1));
nBlocks = size(clustInfo.clusterBlockLocations,2);

% load first cell finding file to get estimate for cell counts and block
% size
cellfile = load(nS.cellFileNameFcn(1,1),'rois');

if isempty(whichClusters)
    whichClusters = 1:nClusts;
end
if isempty(whichBlocks)
    whichBlocks = 1:nBlocks;
end

% load registered peaks and the path to the stack used 
if isempty(refLocSumm)
    refLocSumm = load(nS.referenceLocalizationFileName,'xyzrcoPeak','stackPath');
end

% load stack and get its dimensions
stackFileName = fullfile(PATH_DATA, refLocSumm.stackPath);
stackIs = imageSeries(stackFileName);
stack = squeeze(permute(stackIs.images,[1 2 4 3]));
stackSz = size(stack);

% initialize output variables
xyzrcoClusterPeaks = nan(6, nClusts, nBlocks);
rois(nClusts, nBlocks) = struct('w',[],'x',[],'y',[],'z',[]);

for cc = whichClusters

    fprintf('\nworking on cluster %d, block...',cc);
    
    % list frames associated with this cluster
    thisClustFrames = find(matrixDownsampleSum(clusterFile.clusterSpec(:,cc),1000,1) > 0);
    
    % find localization for these frames
    x=squeeze(refLocSumm.xyzrcoPeak(1,:,thisClustFrames));
    y=squeeze(refLocSumm.xyzrcoPeak(2,:,thisClustFrames));
    z=squeeze(refLocSumm.xyzrcoPeak(3,:,thisClustFrames));
    r=squeeze(refLocSumm.xyzrcoPeak(4,:,thisClustFrames));
    
    % fit a curve to the localization for these frames
    fitZ = fit([x(:) y(:)], z(:),p.Results.zFitStyle,'Robust',p.Results.zFitRobust);
    fitR = fit([x(:) y(:)], r(:),p.Results.rFitStyle,'Robust',p.Results.rFitRobust);
    
    for bb = whichBlocks
        if ~mod(bb,6)
            fprintf('  %d...',bb);
        end

        
        % get positioning and average image for block in the cluster
        [cframes, cx, cy, cIm]= getClusterBlockLocalization(subject,theDate,...
            location,cc,bb);
        
        
        
        
        assert(isequal(thisClustFrames,cframes));
        
        cImSz = size(cIm);
        
        % load the rois for this cluster and block
        try
            cellfile  = load(nS.cellFileNameFcn(cc,bb),'rois');
        catch
            warning(['Warning: aborting without finishing. Could not find'...
                'cluster file for cluster %03d, block %03d'],cc,bb);
            return
        end
        assert(isequal(cImSz, [size(cellfile.rois,1), size(cellfile.rois,2)]));
        
        % set up the neighborhood info struct to constrain possible localizations
        nbrhdInf.xCtr = cx;
        nbrhdInf.yCtr = cy;
        nbrhdInf.zCtr = round(fitZ(cx,cy));
        nbrhdInf.rCtr = fitR(cx,cy);
        
        % set up a binary block location matrix that indicates the subset
        % of the cluster block image that can stand the maximal rotation
        biggestAngle = max(abs(nbrhdInf.rCtr + nbrhdInf.rOffToKeep));
        rotatableSection = makeBlockLocations(cImSz(1),cImSz(2),[1 1], 0, biggestAngle);
        
        % localize this cluster block to the stack
        xyzrcoClusterPeaks(:,cc,bb) = localizeBlockInStackNbrhd(rotatableSection,...
            cIm, stack, nbrhdInf);
        
        bestAngle = xyzrcoClusterPeaks(4,cc,bb);
        bestXCtr = xyzrcoClusterPeaks(1,cc,bb);
        bestYCtr = xyzrcoClusterPeaks(2,cc,bb);
        
        % figure out how big cIm is going to be when rotated appropriately
        [rotH, rotW] = dimAfterRotation(cImSz(1), cImSz(2), bestAngle);
        padSz = ceil([rotH-cImSz(1)  rotW-cImSz(2)]./2);
        
        % pad roi weights with nans so they can rotate
        roiWPadded = padarray(cellfile.rois,padSz,nan,'both');
        roiWLoc = true(size(roiWPadded,1),size(roiWPadded,2));
        %nonNanLoc(~isnan(roiWPadded(:,:,1))) = true;
        
        roiWPaddedSz = size(roiWPadded);
        
        rotatedRoiW = nan(roiWPaddedSz);
        for ww = 1:roiWPaddedSz(3)
            rotatedRoiW(:,:,ww) = rotateAndSelectBlock(roiWPadded(:,:,ww),...
                roiWLoc,bestAngle);
        end
        rois(cc,bb).w = rotatedRoiW;
        rois(cc,bb).x = repmat(bestXCtr + meanCenter(1:roiWPaddedSz(2)),...
            roiWPaddedSz(1),1);
        rois(cc,bb).y = repmat(bestYCtr + meanCenter(1:roiWPaddedSz(1))'...
            , 1,roiWPaddedSz(2));
        rois(cc,bb).z = fitZ(rois(cc,bb).x,rois(cc,bb).y);
    end
end



