% NOT CURRENT
% plots vector field of block movements from movie to ref


load /Volumes/tank/jtb3/Data/j115/2015-12-06__L01__AVERAGE_corrs/J115_2015-12-06_frame001.mat

%% get center of each corr mat
for bb=1:100
    [~,peakZ] = max(best_corr_z(bb,1:50));
    peakCorr = max(reshape(squeeze(C(bb,:,:,peakZ)),[],1)); 
    [ii,jj]=find(squeeze(C(bb,:,:,peakZ))==peakCorr);
    peakLoc(bb,:) = [ii+blk_cntrs.y(bb) jj+blk_cntrs.x(bb)];
    allPeakZ(bb)=peakZ;
end

%% plot
figure(10);clf
cols = jet(52);
for bb=1:100
    %plot(blk_cntrs.x(bb),peakLoc(bb,2),'.-')
    %hold on
    %plot(blk_cntrs.y(bb),peakLoc(bb,1),'.-')
    plot([blk_cntrs.x(bb) peakLoc(bb,2)-100],[blk_cntrs.y(bb) peakLoc(bb,1)-100],'-w')
    hold on
    %plot([blk_cntrs.x(bb) ],[blk_cntrs.y(bb) ],'ok')
    text([blk_cntrs.x(bb) ],[blk_cntrs.y(bb) ],num2str(allPeakZ(bb)),'fontsize',20,'color',cols(allPeakZ(bb),:))
end
    
set(gca,'color',[0 0 0])
axis ij

%% 

figure(11);clf

for bb=1:100
    subplot(1,2,1)
    plot(blk_cntrs.x(bb),peakLoc(bb,2)-100,'ko')
    hold on
    subplot(1,2,2)
    
    plot(blk_cntrs.y(bb),peakLoc(bb,1)-100,'ko')
    hold on
end


%%
bb=52;
[~,peakZ] = max(best_corr_z(bb,:));
figure(13);clf;imagesc(squeeze(C(bb,:,:,peakZ)))
peakCorr = max(reshape(squeeze(C(bb,:,:,peakZ)),[],1)); 
[ii,jj]=find(squeeze(C(bb,:,:,peakZ))==peakCorr);
hold on
plot(jj,ii,'ro')
axis xy
colormap(redbluecmap)
set(gca,'clim',[-.5 .5])


