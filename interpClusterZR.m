function [xyzrcoClusterPeaks cb] = interpClusterZR(subject, theDate, location, varargin)
p = inputParser;
addParamValue(p,'whichClusters',[]);
addParamValue(p,'whichBlocks',[]);
addParamValue(p,'xMargin',15);
addParamValue(p,'yMargin',15);
parse(p,varargin{:});

nbrhdInf.xMargin = p.Results.xMargin;
nbrhdInf.yMargin = p.Results.yMargin;

whichBlocks = p.Results.whichBlocks;
whichClusters = p.Results.whichClusters;

nS = getNameStruct(subject, theDate, location);

% load registered peaks
summ = load(nS.referenceLocalizationFileName,'xyzrcoPeak','stackPath');

% load movie
movieFileName = [nS.avgFinalPrefix '.tif'];
movieInfo = imfinfo(movieFileName);
nFrames = length(movieInfo);
movie = nan(movieInfo(1).Height,movieInfo(1).Width,nFrames);

% get number of blocks and clusters
nBlocks = size(summ.xyzrcoPeak,2);
clusterFile = load(nS.clusterFileName,'clusterSpec');
nClusts = size(clusterFile.clusterSpec,2);

if isempty(whichBlocks)
    whichBlocks=1:nBlocks;
end
if isempty(whichClusters)
    whichClusters = 1:nClusts;
end


nbrhdInf.zOffToKeep = -2:2;
nbrhdInf.rOffToKeep = -1:.25:1;

nbrhdInf.nXYToKeep = 400;

clusterFrames = cell(nClusts,1);

stackFileName = fullfile(PATH_DATA, summ.stackPath);
stackIs = imageSeries(stackFileName);
stack = squeeze(permute(stackIs.images,[1 2 4 3]));

xyzrcoClusterPeaks = nan(6,length(whichClusters),length(whichBlocks));

for cc = whichClusters
    fprintf('\nworking on cluster %d, block...',cc);
    
    % determine which frames are associated with this cluster
    clusterFrames{cc} = find(matrixDownsampleSum(clusterFile.clusterSpec(:,cc),1000,1) > 0);
    x=squeeze(summ.xyzrcoPeak(1,:,clusterFrames{cc}));
    y=squeeze(summ.xyzrcoPeak(2,:,clusterFrames{cc}));
    z=squeeze(summ.xyzrcoPeak(3,:,clusterFrames{cc}));
    r=squeeze(summ.xyzrcoPeak(4,:,clusterFrames{cc}));
    
    fitZ = fit([x(:) y(:)], z(:),'poly11','Robust','Bisquare');
    %fig = figure(11); clf(fig);
    %plot(fZ,[x(:) y(:)], z(:),'parent',subplot(1,2,1,'parent',fig))
    
    fitR = fit([x(:) y(:)], r(:),'poly11','Robust','Bisquare');
    %plot(fR,[x(:) y(:)], r(:),'parent',subplot(1,2,2,'parent',fig))
    
    % fprintf('\nloading movie frames...');
    % movieIs = imageSeries(movieFileName,'whichImages',clusterFrames{cc})

    %meanMovieFrame = squeeze(mean(movieIs.images,4));
    
    %drawnow
    for bb = whichBlocks
        
        
        % let the user know how things are going
        if ~mod(bb,6)
            fprintf('%i... ',bb)
        end
        
        % get x and y position of this cluster
        %%%%%%%%%%%
        % WHAT WE REALLY WANT HERE IS THE WHOLE AV MOVIE FRAME THAT CLUSTER
        % SO THAT WE CAN ALLOW FOR ROTATIONS
        %%%%%%%%%%
        [cf, cx, cy, cIm]= getClusterBlockLocalization(subject,theDate,...
            location,cc,bb);
        assert(isequal(cf,clusterFrames{cc}));  
        
        % set up the neighborhood info struct to constrain possible localizations
        nbrhdInf.xCtr = cx;
        nbrhdInf.yCtr = cy;
        nbrhdInf.zCtr = round(fitZ(cx,cy));
        nbrhdInf.rCtr = fitR(cx,cy);
        
        % set up block location so that rotateAndReselectBlock works
        % properly and create a matched padded version of cIm
        blockLocPad = 100;
        xyBorder = 5;
        thisBlockLoc = false(size(cIm)+[1 1].*blockLocPad*2);
        
        thisBlockLoc((blockLocPad+xyBorder+1):(size(thisBlockLoc,1)-blockLocPad-xyBorder), ...
            (blockLocPad+xyBorder+1):(size(thisBlockLoc,2)-blockLocPad-xyBorder)) = true;
        bInf = getBlockInf(thisBlockLoc);
        cImPadded = zeros(size(thisBlockLoc));
        cImPadded((blockLocPad+1):(size(thisBlockLoc,1)-blockLocPad), ...
            (blockLocPad+1):(size(thisBlockLoc,2)-blockLocPad)) = cIm;
        
        % repeat the process of localizing the block in the stack on the
        % average of the frames in the cluster
        [xyzrcoClusterPeaks(:,cc,bb), cb] = localizeBlockInStackNbrhd(thisBlockLoc,...
            cImPadded, stack, nbrhdInf);
        

    end
    
end

%%


