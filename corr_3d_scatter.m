%load('/Volumes/tank/jtb3/Data/J115/2015-09-17__L01__AVERAGE_corrs/archive/J115_2015-09-17_frame1.mat','C');

subj  = 'J115'
frame = 1;
day   = '2015-09-17'; 

zoffs = -15:15;
offs    = -30:30;
zs = 1:51;
%offix   = find(ismember(off_range,offs));
blockno = 15;

Cblk        = squeeze(C(blockno,:,:,zs));
[~, maxix] = max(Cblk(:));

[maxy maxx maxz] = ind2sub(size(Cblk),maxix);
zs = maxz+zoffs
zs = zs(zs > 0 & zs <= size(C,4));
xix = maxx+offs;
yix = maxy+offs;
Cblk = Cblk(yix,xix,zs);
Cblk_pos                = Cblk;
Cblk_pos(Cblk_pos<0)    = .0010;

[x,y,z] = meshgrid(xix,yix,zs);

%scatter3(x(:),y(:),z(:),70*Cblk_pos(:),Cblk(:),'filled')
slice(Cblk,31,31,find(zs==maxz)); shading flat; colormap redbluecmap(11)
axis image; 
grid off;
colormap redbluecmap;
colorbar;
xlabel('pixels from block center x-coord   ')
ylabel('pixels from block center y-coord   ')
zlabel('stack frame # ')
set(gca, 'XTick',1:10:length(xix),'XTickLabel',off_range(maxx)+offs(1:10:length(xix)),...
    'ZTick',1:5:length(zs),'ZTickLabel',zs(1:5:length(zs)),...
    'YTick',1:10:length(yix),'YTickLabel',off_range(maxy)+offs(1:10:length(yix)));
title(sprintf('%s %s block %i',subj,day,blockno));
keyboard;




%
%
%
%>> Cblk = squeeze(C(15,:,:,:))
%>> figure;imagesc(squeeze(Cblk(10,:,:)))
%>> figure;imagesc(sum(Cblk,3))
%>> figure;imagesc(squeeze(Cblk(120,:,:)))
%>> figure;imagesc(squeeze(Cblk(120,:,:))')
%>> axis image
%>> figure;imagesc(squeeze(Cblk(:,97,:))')
%>> plotImagesSequentially(imageSeries(permute(Cblk,[2 3 1]))
% plotImagesSequentially(imageSeries(permute(Cblk,[2 3 1]))
%                                                           |
%                                                           Error: Expression or statement is incorrect--possibly unbalanced (, {, or [.
%
%                                                           >> plotImagesSequentially(imageSeries(permute(Cblk,[2 3 1])))
%                                                           Undefined function 'imageSeries' for input arguments of type 'double'.
%
%                                                           >> plotImagesSequentially(imageSeries(permute(Cblk,[2 3 1])))
%                                                           >> addpath('/Volumes/jlgauthi/code
%                                                           code backup/       code deletable/    code_overwritten/
%                                                           code copy/         code/
%                                                           >> addpath(genpath('/Volumes/jlgauthi/code/'))
%                                                           >> plotImagesSequentially(imageSeries(permute(Cblk,[2 3 1])))
%                                                           >> plotImagesSequentially(imageSeries(permute(Cblk,[3 2 1])))
%                                                           >> figure;imagesc(max(Cblk,[],3))
