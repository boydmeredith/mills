function plotRegisteredRoisInSection(regRoi,clusterInd,cellDuplicateLists,...
    xyzrcoClusterPeaks,blockInd, stack)

ncolors = 20;
roiThresh = .005;

for xSlice = 80:1:200
    
    title(sprintf('y = %i',xSlice));
    %xlim([100 200])
    %ylim([25 41])
    sliceX = (xSlice:xSlice+3)+100;
    
    
    %     zInRef = xyzrcoClusterPeaks(3,bb,cc);
    
    
    hold on
    roiSetInd = find(cellfun(@(x) length(x)>0, cellDuplicateLists));
    ncolors = 15;
    cmap = hsv(ncolors);
    
    xslice = padarray(squeeze(max(stack(sliceX-100,:,:)))',[0 100],nan,'both');
    yslice = padarray(squeeze(max(permute(stack(:,sliceX-100,:),[2 3 1]))),[0 100],nan,'both');
    imagesc(xslice);
    colormap(bone);
    freezeColors;
    axis image
    
    for mm = 1:length(roiSetInd)
        thisRoiSetInd = roiSetInd(mm);
        roiInd = cellDuplicateLists{thisRoiSetInd};
        
        
        thisRoiColor = cmap(mod(mm,ncolors)+1,:);
        
        
        
        for rr = 1:length(roiInd)
            thisBlockNum = blockInd(roiInd(rr));
            thisClusterNum = clusterInd(roiInd(rr));
            thisClusterZ = xyzrcoClusterPeaks(3,thisClusterNum,thisBlockNum);
            yRoiProfile = find(nansum(regRoi(:,sliceX,roiInd(rr)),2));
            xRoiProfile = find(nansum(regRoi(sliceX,:,roiInd(rr))));
            %keyboard
            thisRoiProfile = xRoiProfile;
            plot(thisRoiProfile,(thisClusterZ+.1*thisClusterNum).*ones(1,length(thisRoiProfile)),'color',thisRoiColor);
            axis ij
        end
        
        %     cf  = load(nS.cellFileNameFcn(cc,bb),'rois');
        %     cSize = size(cf.rois);
        % %     nRois = cSize(3);
        %     xCtrInRef = xyzrcoClusterPeaks(1,bb,cc);
        %     xIndInRef = xyzrcoClusterPeaks(1,bb,cc)+((-cSize/2+1/2):(cSize/2-1/2));
        %     cmap = hsv(nRois);
        %     colormap(cmap);
        
        %     for rr = 1:nRois
        %
        %
        %         thisRoi = find(sum(regRoi(sliceY,cellDuplicateLists{mm}),1));
        %
        %         plot3(thisRoi,ones(size(thisRoi)).*sliceY, ones(size(thisRoi)).*zInRef,'r','linewidth',20);
        %         pause()
        %     end
    end
    pause()
end