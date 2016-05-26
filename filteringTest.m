
xVal = reshape(permute(xyzrcPeak(1,:,1),[2 1]),10,10);
yVal = reshape(permute(xyzrcPeak(2,:,1),[2 1]),10,10);
zVal = reshape(permute(xyzrcPeak(3,:,1),[2 1]),10,10);
rVal = reshape(permute(xyzrcPeak(4,:,1),[2 1]),10,10);
xFilt = medfilt2(xVal,[4 4],'symmetric');
yFilt = medfilt2(yVal,[4 4],'symmetric');
zFilt = medfilt2(zVal,[4 4],'symmetric');
rFilt = medfilt2(rVal,[4 4],'symmetric');

xFilt = xFilt + median(xVal(:)-xFilt(:));
yFilt = yFilt + median(yVal(:)-yFilt(:));
%
figure(21);clf;
subplot(2,2,1);
imagesc([xVal xFilt]);title('x')
subplot(2,2,2);
imagesc([yVal yFilt]);title('y')
subplot(2,2,3);
imagesc([zVal zFilt]);title('z')
subplot(2,2,4);
imagesc([rVal rFilt]);title('r')
colormap(colormapRedBlue)

figure(22);clf
subplot(2,2,1);
plot(xVal(:),xFilt(:),'o');hold on;plot(xlim,xlim,'r-')
subplot(2,2,2);
plot(yVal(:),yFilt(:),'o');hold on;plot(xlim,xlim,'r-')
subplot(2,2,3);
plot(zVal(:),zFilt(:),'o');hold on;plot(xlim,xlim,'r-')
subplot(2,2,4);
plot(rVal(:),rFilt(:),'o');hold on;plot(xlim,xlim,'r-')


%%
ptCloud = pointCloud(xyzrcPeak(1:3,:,1)');
[model,inlierInd,outlierInd] = pcfitplane(ptCloud, 3);
pcshow(select(ptCloud,inlierInd),'markersize',700); axis square
hold all; pcshow(select(ptCloud,outlierInd),'markersize',50);
set(gca,'zdir','reverse','xdir','reverse','ydir','reverse')
xlabel('x'); ylabel('y'); zlabel('z')