function [clustCtrX, clustCtrY, bestZ, bestR] = interpClusterZR(subject, theDate, location, varargin)
p = inputParser;
addOptional(p,'whichClusters',[]);
addOptional(p,'whichBlocks',[]);
parse(p,varargin{:});

whichBlocks = p.Results.whichBlocks;
whichClusters = p.Results.whichClusters;

nS = getNameStruct(subject, theDate, location);

% load registered peaks
summ = load(nS.referenceLocalizationFileName,'xyzrcoPeak');

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

% initialize output variables
bestZ = nan(nClusts,nBlocks); bestR = nan(nClusts,nBlocks);
clustCtrX = nan(nClusts,nBlocks);
clustCtrY = nan(nClusts,nBlocks);

for cc = whichClusters
    fprintf('\nworking on cluster %d, block...',cc);
    
    % determine which frames are associated with this cluster
    clusterFrames{cc} = find(matrixDownsampleSum(clusterFile.clusterSpec(:,cc),1000,1) > 0);
    x=squeeze(summ.xyzrcoPeak(1,:,clusterFrames{cc}));
    y=squeeze(summ.xyzrcoPeak(2,:,clusterFrames{cc}));
    z=squeeze(summ.xyzrcoPeak(3,:,clusterFrames{cc}));
    r=squeeze(summ.xyzrcoPeak(4,:,clusterFrames{cc}));
    
    fZ = fit([x(:) y(:)], z(:),'poly11','Robust','Bisquare');
    %fig = figure(11); clf(fig);
    %plot(fZ,[x(:) y(:)], z(:),'parent',subplot(1,2,1,'parent',fig))
    
    fR = fit([x(:) y(:)], r(:),'poly11','Robust','Bisquare');
    %plot(fR,[x(:) y(:)], r(:),'parent',subplot(1,2,2,'parent',fig))
    
    %drawnow
    for bb = whichBlocks
        
        % let the user know how things are going
        if ~mod(bb,6)
            fprintf('%i... ',bb)
        end
        
        % get x and y position of this cluster
        [cf, clustCtrX(cc,bb), clustCtrY(cc,bb) ]= getClusterBlockLocalization(subject,theDate,location,cc,bb);
        assert(isequal(cf,clusterFrames{cc}));

        bestZ(cc,bb) = fZ(clustCtrX(cc,bb), clustCtrY(cc,bb));
        bestR(cc,bb) = fR(clustCtrX(cc,bb), clustCtrY(cc,bb));
        
    end
end
%%


