%%



%% plot correlation stack
figure;
for z=30:51; hold on; thisXC = normxcorr2(block1,stack(:,:,z)); xc(:,:,z) = thisXC(bInf.indY,bInf.indX)       ;surf(xx,yy,2*z.*ones(size(xx)),xc(:,:,z));end; axis image; colormap(redbluecmap); axis off; figure(gcf); shading flat; axis image; set(gca,'zdir','rev'); axis ij
axis ij
set(gcf,'position',thisPos); set(gca,'view',thisView); set(gca,'zdir','rev')
 crange = get(gca,'clim');
 %%
 figure(20); clf
imagesc(squeeze(xc(:,:,z))); axis image; axis ij;
  
  colorbar('southoutside'); colormap redbluecmap
    caxis(crange)
    axis off

 colormap redbluecmap
 figure(21); clf
 imagesc(squeeze(xc(:,y,:))'); axis image; axis ij;
  axis off
  caxis(crange)
  
  %%
   figure(20);
  colorbar('southoutside'); figure(gcf)

  
  %%
  xc30 = normxcorr2(block1,stack(:,:,30));
  [~, maxIx30] = max(xc30(:));
  [y x z] = ind2sub(size(xc30),maxIx30);
  foundXInd = (x-blockWidth+1):x;
  foundYInd = (y-blockHeight+1):y;
  figure(22); clf
  imagesc(stack(foundYInd,foundXInd,30));
  axis ij
  axis image
  colormap bone
 
  
  %%
  figure(30); clf
  ztolookat=30:34;
  for zz = 1:length( ztolookat)
      thisZ = ztolookat(zz);
      subplot(length(ztolookat),1,zz);
      
      imagesc(stack(foundYInd,foundXInd,thisZ));
      axis image
      colormap bone
  end
  
  %
 
%% plot stack
 figure; 
 for z=1:51; hold on;      ;surf(xxS,yyS,2*z.*ones(size(xxS)),normalizeToZeroOne(double(stack(:,:,z))));end; axis image; colormap(bone); axis off; figure(gcf); shading flat; axis image; set(gca,'zdir','rev'); axis ij
 axis ij
 set(gcf,'position',thisPos); set(gca,'view',thisView); set(gca,'zdir','rev')
 %% plot slice
 figure; 
 slice(xc,exSumm1103.xyzrcoPeak(1,1,1)-89/2,exSumm1103.xyzrcoPeak(2,1,1)-90/2, exSumm1103.xyzrcoPeak(3,1,1))
 axis ij
 set(gcf,'position',thisPos); set(gca,'view',thisView); set(gca,'zdir','rev')
 
 
 
 