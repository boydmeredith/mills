
% load
jlgDataDir = '/Volumes/tank/jlgauthi/Data';
subj = 'J115';
subjDir = fullfile(jlgDataDir,subj);
stackPath = fullfile(subjDir,'2015-09-25__post_stack_002_AVERAGE.gif');
moviePath = fullfile(subjDir,'2015-12-06__L01__AVERAGE.gif');
movieFrame = imread(moviePath,1);
%%

[movieFrameHeight movieFrameWidth] = size(movieFrame);
blockLocations = makeBlockLocations(movieFrameHeight, movieFrameWidth, 10, .2, 5);


%%
thisBlockLoc = blockLocations(:,:,1);
[bb] = getBlockInf(thisBlockLoc);

bb.height = length(bbY);
bb.width = length(bbX);

bb.ctrY = mean(bbY);
bb.ctrX = mean(bbX);

rotationAngles = -10:1:10;
nRotAngles = length(rotationAngles);



blockLocMarked = thisBlockLoc;
blockLocMarked(:,[floor(bb.ctrX) ceil(bb.ctrX)]) = false;
blockLocMarked([floor(bb.ctrY) ceil(bb.ctrY)],:) = false;
figure(1); clf; imagesc(blockLocMarked);

imTF = zeros(bb.height,bb.width,nRotAngles);

blockToRotate = movieFrame.*uint8(blockLocMarked);%blockLocMarked;

for rr = 1:nRotAngles
    r=rotationAngles(rr);
    T = maketform('affine',[cosd(r) sind(r) 0; -sind(r) cosd(r) 0; 0 0 1]);
    imTF(:,:,rr)=(imtransform(blockToRotate,T,'bilinear',...
        'udata',[1 movieFrameWidth]  - bb.ctrX,...
        'vdata',[1 movieFrameHeight] - bb.ctrY,...
        'xdata',([1 bb.width]) -(bb.width+1)/2,...
        'ydata',([1 bb.height])-(bb.height+1)/2,...
        'xyscale',1));
    
end

figure(2);colormap(colormapRedBlue)
plotImagesSequentially(imageSeries(imTF));
figure(3);
imagesc(sum(imTF,3)); colormap(colormapRedBlue)