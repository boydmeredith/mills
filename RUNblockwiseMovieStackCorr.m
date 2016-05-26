%% load
jlgDataDir = '/Volumes/tank/jlgauthi/Data';
subj = 'J115';
subjDir = fullfile(jlgDataDir,subj);
stackName = 'reference_stack_2015-09-25.tif';
stackPath = fullfile(subjDir,stackName);
movieName = '2015-12-06__L01__AVERAGE.tif';
moviePath = fullfile(subjDir,movieName);
frameOneCorrsPath = fullfile(subjDir,'2015-12-06__L01__AVERAGE_corrs/J115_2015-12-06_frame001.mat');

%%
stackInf = imfinfo(stackPath);
stack(stackInf(1).Height,stackInf(1).Width,length(stackInf)) = 0;
for ss = 1:length(imfinfo(stackPath))
    stack(:,:,ss) = imread(stackPath,ss);
end
stack = cropStack(stack);
[sY, sX, sZ] = size(stack);

movieInf    = imfinfo(moviePath);
movieHeight = movieInf(1).Height;
movieWidth  = movieInf(1).Width;
movieLength = length(movieInf); 
movie(movieHeight,movieWidth,movieLength) = 0;
for mm = 1:length(imfinfo(moviePath))
    movie(:,:,mm) = imread(moviePath,mm);
end
%% run blockwise correlations fn
blockwiseMovieStackCorr(subj,movieName,stackName,'loadedMovie',movie,'loadedStack',stack,...
    'whichSlices',[],'coarseRotWindowRange',8, 'showFigs','on');