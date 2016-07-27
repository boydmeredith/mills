%% add code to path
addJeffCode2Path
%%
subj = 'J114';

movieDate = '2015-09-18';
%%
referenceDateSt = struct('J114','2015-11-30',...
                         'J115', '2015-09-25',...
                         'J116', '2015-09-28',...
                         'J117', '2015-09-25',...
                         'J118', '2015-10-01',...
                         'J122', '2015-09-27',...
                         'J123', '2015-09-25');
%
if strcmp(movieDate,'stackDate')
    movieDate = referenceDateSt.(subj);
end

jlgDataDir = '/Volumes/tank/jlgauthi/Data';
if ~exist(jlgDataDir)
    jlgDataDir = '/jukebox/tank/jlgauthi/Data';
end
subjDir = fullfile(jlgDataDir,subj);
stackName = ['reference_stack_' referenceDateSt.(subj) '.tif'];

stackPath = fullfile(subjDir,stackName);

movieName = sprintf('%s__L01__AVERAGE.tif',movieDate);
moviePath = fullfile(subjDir,movieName);


stackInf = imfinfo(stackPath);
stack = zeros(stackInf(1).Height,stackInf(1).Width,length(stackInf));
for ss = 1:length(imfinfo(stackPath))
    stack(:,:,ss) = imread(stackPath,ss);
end
stack = cropStack(stack);
[sY, sX, sZ] = size(stack);

movieInf    = imfinfo(moviePath);
movieHeight = movieInf(1).Height;
movieWidth  = movieInf(1).Width;
movieLength = length(movieInf); 
movie = zeros(movieHeight,movieWidth,movieLength);
for mm = 1:movieLength
    movie(:,:,mm) = imread(moviePath,mm);
end
%
%
%%
profile off
profile on
blockwiseMovieStackCorr(subj,movieDate,...
    'loadedMovie',movie,'loadedStack',stack,...
    'whichSlices',[], 'showFigs','on', ...
    'useSavedSearchRange', true, 'useSavedSearchRangeEitherWay', false,...
    'whichBlocks',[],...
    'whichFrames',[],...
       'summarySaveName','');
profile off
profile report 