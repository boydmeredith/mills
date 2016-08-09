
subject='J115';
theDate = '2015-09-25';
location = 'L01';
nS = getNameStruct(subject,theDate,location);
whichClusters = 1;
whichBlocks = 1+ [1 7 13 19 25 ]; %29];

%%
bb = 1;
theDate='2015-11-03'
fprintf('getting peaks')
[xyzrcoClusterPeaks1103 cb] = interpClusterZR(subject, theDate, location, ...
    'whichClusters', whichClusters, 'whichBlocks',whichBlocks)

xyzrcoClusterPeaks1103= xyzrcoClusterPeaks1103(:,whichClusters,:);
%
%plotImagesSequentially(imageSeries(regRoi(100:200,100:200,cellDuplicateLists{15})));
fprintf('registering rois')
[regRoi1103, clusterInd1103, blockInd1103] = registerRois(subject,theDate,location,...
    xyzrcoClusterPeaks1103,stack,'whichClusters', whichClusters, 'whichBlocks',whichBlocks);
fprintf('finding duplicates')
[uniqueCellList1103,cellDuplicateLists1103,extras1103] = findRoiDuplicates(regRoi1103);

%%

bb = 1;
theDate='2015-09-25';
fprintf('getting peaks')
[xyzrcoClusterPeaks0925 cb] = interpClusterZR(subject, theDate, location, ...
    'whichClusters', whichClusters, 'whichBlocks',whichBlocks)

xyzrcoClusterPeaks0925= xyzrcoClusterPeaks0925(:,whichClusters,:);
%
%plotImagesSequentially(imageSeries(regRoi(100:200,100:200,cellDuplicateLists{15})));
fprintf('registering rois')
[regRoi0925, clusterInd0925, blockInd0925] = registerRois(subject,theDate,location,...
    xyzrcoClusterPeaks0925,stack,'whichClusters', whichClusters, 'whichBlocks',whichBlocks);
fprintf('finding duplicates')
[uniqueCellList0925,cellDuplicateLists0925,extras0925] = findRoiDuplicates(regRoi0925);

%%
figure(10); clf
plotRegisteredRoisInSection(regRoi0925,clusterInd0925,cellDuplicateLists0925,...
    xyzrcoClusterPeaks0925,blockInd0925, stack)
%%
figure(11); clf
plotRegisteredRoisInSection(regRoi1103,clusterInd1103,cellDuplicateLists1103,...
    xyzrcoClusterPeaks1103,blockInd1103, stack)
%%
figure(7); clf
plotRegRoiInVol(regRoi0925,clusterInd0925,cellDuplicateLists0925,...
    xyzrcoClusterPeaks0925,blockInd0925, stack)


%%
%%


%%
cf  = load(nS.cellFileNameFcn(1,bb),'rois');
cf2  = load(nS.cellFileNameFcn(2,bb),'rois');
cat(3,cf.rois,cf2.rois);