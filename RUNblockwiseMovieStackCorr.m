%% add code to path
addJeffCode2Path
%%

jlgDataDir = '/Volumes/tank/jlgauthi/Data';
if ~exist(jlgDataDir)
    jlgDataDir = '/jukebox/tank/jlgauthi/Data';
end
referenceDateSt = struct('J114','2015-11-30',...
                         'J115', '2015-09-25',...
                         'J116', '2015-09-28',...
                         'J117', '2015-09-25',...
                         'J118', '2015-10-01',...
                         'J122', '2015-09-27',...
                         'J123', '2015-09-25');
%
subj = 'J118'; 

subjDir = fullfile(jlgDataDir,subj);
stackName = sprintf('reference_stack_%s.tif',referenceDateSt.(subj));
stackPath = fullfile(subjDir,stackName);
movieName = [referenceDateSt.(subj) '__L01__AVERAGE.tif'];
moviePath = fullfile(subjDir,movieName);



%
% stackInf = imfinfo(stackPath);
% stack(stackInf(1).Height,stackInf(1).Width,length(stackInf)) = 0;
% for ss = 1:length(imfinfo(stackPath))
%     stack(:,:,ss) = imread(stackPath,ss);
% end
% stack = cropStack(stack);
% [sY, sX, sZ] = size(stack);
% 
% movieInf    = imfinfo(moviePath);
% movieHeight = movieInf(1).Height;
% movieWidth  = movieInf(1).Width;
% movieLength = length(movieInf); 
% movie(movieHeight,movieWidth,movieLength) = 0;
% for mm = 1:length(imfinfo(moviePath))
%     movie(:,:,mm) = imread(moviePath,mm);
% end
% %
% run blockwise correlations fn
%%
profile on
blockwiseMovieStackCorr(subj,movieName,stackName,...%'loadedMovie',movie,'loadedStack',stack,...
    'whichSlices',[],...%'coarseRotWindowRange',10, 'showFigs','on', ...
    ...%'nbrhdXMargin',10,'nbrhdYMargin',10,'nXYToKeep', 300,...
    'useSavedSearchRange', true, 'useSavedSearchRangeEitherWay', true,...
    'whichBlocks',[],...
    'whichFrames',[70],...
    'dataDir',jlgDataDir);
profile off
