% load
jlgDataDir = '/Volumes/tank/jlgauthi/Data';
subj = 'J115';
subjDir = fullfile(jlgDataDir,subj);
stackPath = fullfile(subjDir,'2015-09-25__post_stack_002_AVERAGE.gif');
moviePath = fullfile(subjDir,'2015-12-06__L01__AVERAGE.gif');
frameOneCorrsPath = fullfile(subjDir,'2015-12-06__L01__AVERAGE_corrs/J115_2015-12-06_frame001.mat');
%%
stack = imread(stackPath,'Frames','all');
stack = cropStack(stack);
[sY, sX, sZ] = size(stack);
movieFrameOne = imread(moviePath,1);
[mY, mX, mZ] = size(movieFrameOne);
%%
stackSlice = stack(:,:,31);

%% divide into blocks
angle = 5;
nBlockSpan = 10;
overlap = .1;
bb = 13;
[blockLocBin] = makeBlockLocations(mY, mX, nBlockSpan, overlap, angle);
[bbY, bbX] = getBlockIx(blockLocBin(:,:,bb));
figure(2); clf
block = movieFrameOne(bbY,bbX); 
imagesc(block);
colormap bone

rotBlock = rotateAndReselectBlock(movieFrameOne, blockLocBin(:,:,bb), 0);
imshowpair(block,rotBlock,'falsecolor')


%% best Z overall has no clear peak
wholeFrameCorr = zeros(1,sZ);
for zz = 1:sZ
    c = normxcorr2(movieFrameOne,stack(:,:,zz));
    wholeFrameCorr(zz) = max(c(:));
end

figure(3); clf
plot(wholeFrameCorr,'.k')



%% best Z overall for specific block
figure(2); clf
[h, hM] = makeSubplots(2, nBlockSpan, nBlockSpan, .01,.01, [0 0 1 1]);

for bb = 1:nBlockSpan^2
    [bbY, bbX] = getBlockIx(blockLocBin, bb);
    smallZ = 1:.1:sZ;
    
    block = movieFrameOne(bbY,bbX);
    maxCorrByZ = zeros(1,sZ);
    for zz = 1:sZ
        c = normxcorr2(block,stack(:,:,zz));
        maxCorrByZ(zz) = max(c(:));
    end
    
    f = fit([1:sZ].',maxCorrByZ.','gauss2');
    [fmax, fmaxix] = max(f(smallZ));
    
    axes(h(bb));
    
    plot([smallZ(fmaxix) smallZ(fmaxix)], [.2 .9], '--', 'color',[.7 .7 .7])
    hold on
    plot(f,1:sZ,maxCorrByZ)
    xlim([1 sZ])
    ylim([.2 .9])
    set(gca, 'xcolor','w','ycolor','w')
    legend('off')
    box off
    drawnow
    
    
end

%% block 24 has a clear peak and plenty of rotation, so we'll use this as a test case
bb = 24;
neighborhood = 20;
[bbY, bbX] = getBlockIx(blockLocBin, bb);
block = movieFrameOne(bbY, bbX);
%%
maxCorrByZ = zeros(1,sZ);
for zz = 1:sZ
    c = normxcorr2(block,stack(:,:,zz));
    maxCorrByZ(zz) = max(c(:));
end

% find
f = fit([1:sZ].',maxCorrByZ.','gauss2');

figure(4); clf
plot([smallZ(fmaxix) smallZ(fmaxix)], [.2 .9], '--', 'color',[.7 .7 .7])
hold on
plot(f,1:sZ,maxCorrByZ)
bestZ = round(smallZ(fmaxix));
%%
bbYS = (bbY(1)-neighborhood):(bbY(end)+neighborhood);
bbXS = (bbX(1)-neighborhood):(bbX(end)+neighborhood);
stackBlockArea = stack(bbYS,bbXS,bestZ);

%%
figure(5); clf
imshowpair(block, stackBlockArea, 'montage')


%% Let's try finding the best z for a block and then looking for the best 
%% rotation within a window in the best z
movieFrameCenter = movieFrameOne(200:300,300:400);
bestCorrByZ = zeros(1,51);
figure(12); clf
for zz = 1:51
    
    
    stackSlice = stack(:,:,zz);
    c = normxcorr2(movieFrameCenter,stackSlice);
    bestCorrByZ(zz) = max(c(:));
    [peakY peakX] = ind2sub(size(c), find(c == max(c(:))));
    offY = peakY-size(movieFrameCenter,1)+1;
    offX = peakX-size(movieFrameCenter,2)+1;
    paddedFrame = padarray(movieFrameCenter,[max(0,offY) max(0,offX)],0);
    paddedStackSlice = padarray(stackSlice,[max(0,-offY) max(0,-offX)],0);
    
    subplot(2,2,1)
    imshowpair(paddedFrame,paddedStackSlice,'montage')
    subplot(2,2,2)
    imshowpair(paddedFrame,paddedStackSlice,'falsecolor')
    subplot(2,2,3:4)
    scatter(zz, bestCorrByZ(zz),'ko')
    xlabel('z')
    hold on
    drawnow
end
[bestCorr bestZ] = max(bestCorrByZ);
stackSlice = stack(:,:,bestZ);

c = normxcorr2(movieFrameCenter,stackSlice);
bestCorr = max(c(:));
[peakY peakX] = ind2sub(size(c), find(c == max(c(:))));
offY = peakY-size(movieFrameCenter,1)+1;
offX = peakX-size(movieFrameCenter,2)+1;
paddedFrame = padarray(movieFrameCenter,[max(0,offY) max(0,offX)],0);
paddedStackSlice = padarray(stackSlice,[max(0,-offY) max(0,-offX)],0);
subplot(2,2,1)
imshowpair(paddedFrame,paddedStackSlice,'montage')
subplot(2,2,2)
imshowpair(paddedFrame,paddedStackSlice,'falsecolor')
subplot(2,2,3:4)
scatter(bestZ, bestCorrByZ(bestZ),'k.')
drawnow

%%
figure(12);clf
stackSliceRegion = stackSlice(150:350,250:450);
bestCorrByR = zeros(size(rotations));
for rr = 1:length(rotations)
    movieFrameRot = imrotate(movieFrameOne, rotations(rr));
    movieFrameCenter = movieFrameRot(200:300,300:400);

    c = normxcorr2(movieFrameCenter,stackSliceRegion);
    bestCorrByR(rr) = max(c(:));
    
    [peakY peakX] = ind2sub(size(c), find(c == max(c(:))));
    offY = peakY-size(movieFrameCenter,1)+1;
    offX = peakX-size(movieFrameCenter,2)+1;
    paddedFrame = padarray(movieFrameCenter,[max(0,offY) max(0,offX)],0);
    paddedStackSlice = padarray(stackSliceRegion,[max(0,-offY) max(0,-offX)],0);
    subplot(2,2,1)
    imshowpair(paddedFrame,paddedStackSlice,'montage')
    subplot(2,2,2)
    imshowpair(paddedFrame,paddedStackSlice,'falsecolor')
    drawnow
    subplot(2,2,3:4); hold on
    scatter(rotations(rr),bestCorrByR(rr),'ko');
    xlabel('angle')
    
end

[bestRotCorr bestRotIx] = max(bestCorrByR);
bestRot = rotations(bestRotIx);

movieFrameZRot = imrotate(movieFrameOne, bestRot);
% it's not good to grab a different section of the image for each
% correlation and then try to compare those correlations.

%% try using detectSURFfeatures -- BAD
ptsOriginal = detectSURFFeatures(stackBlockArea);
ptsDistorted = detectSURFFeatures(block);

[featuresOriginal,   validPtsOriginal]  = extractFeatures(stackBlockArea,  ptsOriginal);
[featuresDistorted, validPtsDistorted]  = extractFeatures(block, ptsDistorted);

indexPairs = matchFeatures(featuresOriginal, featuresDistorted);

matchedOriginal  = validPtsOriginal(indexPairs(:,1));
matchedDistorted = validPtsDistorted(indexPairs(:,2));

figure;
showMatchedFeatures(stackBlockArea,block,matchedOriginal,matchedDistorted);
title('Putatively matched points (including outliers)');
%% Using imregcorr to estimate translation AND rotation within neighborhood -- BAD
block = imgaussfilt(block,.25);
stackBlockArea = imgaussfilt(stackBlockArea,.25);
tformEstimate = imregcorr(block,stackBlockArea,'rigid');
Rfixed = imref2d(size(stackBlockArea));
blockReg = imwarp(block,tformEstimate,'OutputView',Rfixed);

figure(2); clf
subplot(3,1,1); imshowpair(stackBlockArea,block,'montage');
subplot(3,1,2); imshowpair(stackBlockArea,blockReg,'montage');
subplot(3,1,3); imshowpair(stackBlockArea,blockReg,'falsecolor');
%% Try using imregcorr with whole image - also BAD
stackSlice = stack(:,:,bestZ);
stackSliceDown = imgaussfilt(stackSlice,2);
movieFrameOneDown = imgaussfilt(movieFrameOne,2);
tformEstimate = imregcorr(movieFrameOneDown,stackSliceDown,'rigid');

Rfixed = imref2d(size(stackSliceDown));


frameReg = imwarp(movieFrameOneDown,tformEstimate,'OutputView',Rfixed);

figure(2); clf
subplot(3,1,1); imshowpair(stackSliceDown,movieFrameOneDown,'montage');
subplot(3,1,2); imshowpair(stackSliceDown,frameReg,'montage');
subplot(3,1,3); imshowpair(stackSliceDown,frameReg,'falsecolor');















