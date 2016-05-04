
% experiment with the method of dividing the image into blocks 

jlgDataDir = '/Volumes/tank/jlgauthi/Data/';
exampleMovie = fullfile(jlgDataDir, 'J115/2015-09-25/L01/ac_001_073.mat');
m = load(exampleMovie);

% get safe zone of motioned corrected 1000 frame mean
safeMean = findSafeZone(m.shifts, m.finalMean);

% TODO: what if the image is actually a 3000 frame mean???


%%
% divide safe mean into blocks
nBars = 10;
percentOverlap = .3;
[safeH, safeW] = size(safeMean); 
[startsH, endsH, lengthH] = divideIntoBlocks(safeH, nBars, percentOverlap);
[startsW, endsW, lengthW] = divideIntoBlocks(safeW, nBars, percentOverlap);

    
% record binary matrix of block locations
[~,hM] = makeSubplots(1,nBars,nBars,.01,.01,[0 0 1 1]);
blockLocations = makeBlockLocations(safeMean, nBars, percentOverlap);
bb = 1;
for wb = 1:nBars
    for hb = 1:nBars
        block = (safeMean(logical(blockLocations(:,:,bb))));
        height = sum(max(blockLocations(:,:,bb),[],1));
        width = sum(max(blockLocations(:,:,bb),[],2));
        axes(hM(hb,wb));
        imagesc(reshape(block,[height width]));
        set(gca,'XTick',[]);
        set(gca,'YTick',[]);
        bb = bb+1;
    end
end
figure(2); colormap bone
imagesc(sum(blockLocations,3))
