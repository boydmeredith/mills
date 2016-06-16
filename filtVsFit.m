
[xx, yy]  = meshgrid(1:10,1:10);
figure(10); clf
set(gcf, 'position', [10  1   1611  947])

nRSTD = 8;

subplot(3,3,1)
[fX, gofX, outX]  = fit([xx(:), yy(:)], xyzrPrior(1,:)', 'poly11','Robust','Bisquare')
%fXLin = fit([xx(:), yy(:)], xyzrPrior(1,:)', 'poly11','Robust','LAR')
plot(fX,[xx(:),yy(:)],xyzrPrior(1,:)')
set(gca,'ydir', 'reverse', 'view', [-13.9000    7.6000]);
colormap(colormapRedBlue)
axis vis3d
xlabel('x block #'); ylabel('y block #'); zlabel('x center (pixels)');
subplot(3,3,4)
% plot(fXLin,[xx(:),yy(:)],xyzrPrior(1,:)')
% set(gca,'ydir', 'reverse', 'view', [-13.9000    7.6000]);
% axis vis3d
% xlabel('x block #'); ylabel('y block #'); zlabel('x center (pixels)');
hist(outX.residuals,5000)
hold on
plot([std(outX.residuals) std(outX.residuals)],get(gca,'ylim'),'r',...
    -[std(outX.residuals) std(outX.residuals)],get(gca,'ylim'),'r',...
    nRSTD*[robustSTD(outX.residuals) robustSTD(outX.residuals)],get(gca,'ylim'),'m',...
    -nRSTD*[robustSTD(outX.residuals) robustSTD(outX.residuals)],get(gca,'ylim'),'m')
% imagesc(abs(reshape(outX.residuals,10,10)) > 3.5*robustSTD(outX.residuals)); title('residuals > robust std')
% axis image; colorbar
subplot(3,3,7)

outliersX = abs(outX.residuals) > nRSTD*robustSTD(outX.residuals);
imagesc(abs(reshape(outX.residuals,10,10)) > nRSTD*robustSTD(outX.residuals)); title('residuals > std')
axis image; colorbar




subplot(3,3,2)
[fY, gofY, outY] = fit([xx(:), yy(:)], xyzrPrior(2,:)', 'poly11','Robust','Bisquare')
%fYLin = fit([xx(:), yy(:)], xyzrPrior(2,:)', 'poly11','Robust','LAR')
plot(fY,[xx(:),yy(:)],xyzrPrior(2,:)')
set(gca,'ydir', 'reverse','view',[-84.8000   10.4000]);
axis vis3d
xlabel('x block #'); ylabel('y block #'); zlabel('y center (pixels)');
subplot(3,3,5)
% plot(fYLin,[xx(:),yy(:)],xyzrPrior(2,:)')
% set(gca,'ydir', 'reverse','view',[-84.8000   10.4000]);
% axis vis3d
% xlabel('x block #'); ylabel('y block #'); zlabel('y center (pixels)');    
hist(outY.residuals,5000)
hold on
plot([std(outY.residuals) std(outY.residuals)],get(gca,'ylim'),'r',...
    -[std(outY.residuals) std(outY.residuals)],get(gca,'ylim'),'r',...
    nRSTD*[robustSTD(outY.residuals) robustSTD(outY.residuals)],get(gca,'ylim'),'m',...
    -nRSTD*[robustSTD(outY.residuals) robustSTD(outY.residuals)],get(gca,'ylim'),'m')
% imagesc(abs(reshape(outY.residuals,10,10)) > 3.5*robustSTD(outY.residuals)); title('residuals > robust std')
% axis image; colorbar
subplot(3,3,8)
imagesc(abs(reshape(outY.residuals,10,10)) > nRSTD*robustSTD(outY.residuals)); title('residuals > std')
axis image; colorbar
outliersY = abs(outY.residuals) > nRSTD*robustSTD(outY.residuals);


outliers = outliersX | outliersY;
%%

subplot(3,3,3)
fZExclude = fit([xx(~outliers), yy(~outliers)], xyzrPrior(3,~outliers)', 'loess','Robust','off')
fZInclude = fit([xx(:), yy(:)], xyzrPrior(3,:)', 'loess','Robust','off')

figure(13); clf

fRExclude = fit([xx(~outliers), yy(~outliers)], xyzrPrior(4,~outliers)', 'poly11','Robust','off')
fRInclude = fit([xx(:), yy(:)], xyzrPrior(4,:)', 'poly11','Robust','off')
%plot(fR,[xx(:),yy(:)],xyzrPrior(4,:)')
colormap(colormapRedBlue)
subplot(2,2,1)
imagesc(fZExclude(xx,yy)-reshape(xyzrPrior(3,:),10,10)); colorbar; caxis([-20 20])
subplot(2,2,3)
imagesc(fZInclude(xx,yy)-reshape(xyzrPrior(3,:),10,10)); colorbar; caxis([-20 20])
subplot(2,2,2)
imagesc(fRExclude(xx,yy)-reshape(xyzrPrior(4,:),10,10)); colorbar; caxis([-6 6])
subplot(2,2,4)
imagesc(fRInclude(xx,yy)-reshape(xyzrPrior(4,:),10,10)); colorbar; caxis([-6 6])

%%
ptCloud = pointCloud([xx(:), yy(:), xyzrPrior(1,:)']);
[model,inlierInd,outlierInd] = pcfitplane(ptCloud, 100);

figure(12); clf
pcshow(select(ptCloud,inlierInd),'markersize',5000); 
hold on
axis square
pcshow(select(ptCloud,outlierInd),'markersize',5000);



xlabel('x'); ylabel('y'); zlabel('z')   

set(gca,'zdir','reverse','xdir','reverse','ydir','normal')
axis square
grid off
shading flat


plot(model);