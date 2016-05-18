
%% QUESTIONS FOR JEFF: 
% find out if I have shifts right (i.e., x is element 2)
% what if the image is actually a 3000 frame mean???
% find out what has been done to the 

%% load example image
jlgDataDir = '/Volumes/tank/jlgauthi/Data/';
movPath = fullfile(jlgDataDir, 'J115/2015-09-25/L01/ac_001_073.mat');
mov = load(movPath);
% get safe zone of motioned corrected 1000 frame mean
meanImg = mov.finalMean; %findSafeZone(mov.shifts, mov.finalMean);

%% load corresponding stack
stackPath = fullfile(jlgDataDir, 'J115/2015-09-25__post_stack_002_AVERAGE.tif');     
stack = load_stack(stackPath);

%% divide safe mean into blocks and show the divisions
nBars = 10;
percentOverlap = .3;
[safeH, safeW] = size(meanImg); 
rotAngle = 10;
    
% record binary matrix of block locations
[~,hM] = makeSubplots(1,nBars,nBars,.01,.01,[0 0 1 1]);
blockLocations = makeBlockLocations(safeH,safeW, nBars, percentOverlap, rotAngle);
bb = 1;
for wb = 1:nBars
    for hb = 1:nBars
        [bby bbx] = getBlockIx(squeeze(blockLocations(:,:,bb)));
        block = (meanImg(bby,bbx));
        axes(hM(hb,wb));
        imagesc(block);
        set(gca,'XTick',[]);
        set(gca,'YTick',[]);
        bb = bb+1;
    end
end
figure(2); colormap bone
imagesc(sum(blockLocations,3))

%% experiment with cross correlation and rotation
stackSlice = stack(:,:,30);
blockno = 12;
block = meanImg(find(max(blockLocations(:,:,blockno),[],1)),...
    find(max(blockLocations(:,:,blockno),[],2)));

height = sum(max(blockLocations(:,:,blockno),[],1));
width = sum(max(blockLocations(:,:,blockno),[],2));

corrMat = normxcorr2(block, stackSlice);
[maxCorr maxIx] = max(corrMat(:));
[maxY maxX] = ind2sub(size(corrMat),maxIx);

stackY = maxY - size(block,1) + 1;
stackX = maxX - size(block,2) + 1;
stackBlock = stackSlice([stackY:stackY+height-1],[stackX:stackX+width-1]);

figure(3); clf
subplot(1,4,1);
imagesc(block);
axis image

subplot(1,4,2);
imagesc(stackBlock);
axis image

subplot(1,4,3);
imshowpair(block, stackBlock, 'falsecolor');%, 'ColorChannels','red-cyan');
axis image

subplot(1,4,4);
imagesc(corrMat);
axis image
hold on
%scatter(maxX,maxY,'ro','markersize',5);
colormap bone
caxis([-.5 .8])

%% Using imregcorr to estimate translation AND rotation
tformEstimate = imregcorr(block,stackSlice,'rigid');
Rfixed = imref2d(size(StackSlice));
blockReg = imwarp(block,tformEstimate,'OutputView',Rfixed);

figure(2); subplot(3,1,1); imshowpair(stackSlice,block,'montage');
subplot(3,1,2); imshowpair(stackSlice,blockReg,'montage');
subplot(3,1,3); imshowpair(StackSlice,blockReg,'falsecolor');







