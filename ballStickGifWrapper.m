function hArr = ballStickGifWrapper(subj, movieDate, location, whichFrames, viewPoint, doSave, doShow)
% function hArr = ballStickGifWrapper(subj, movieDate, location, whichFrames, viewPoint, doSave)
% 
% Plot blockwise alignments for the specified frames from a registered
% moviewith each blocks' center shown as a point in x,y,z space, with color
% corresponding to rotation and correlation magnitude shown as point size.
%
% viewPoint - optional argument. Shows bird's eye view if 'birdsEye', shows
% side view otherwise.
%
% output:
% hArr        an array of figures corresponding to each frame of the movie.
%             they will not be visible until the user asks for them
%
% Saves a gif if desired



if isempty(location), location = 'L01'; end
if isempty(doSave), doSave = false; end

theRefLocDir = referenceLocalizationDir(subj, movieDate, location);
s = load(fullfile(theRefLocDir, 'summary.mat'));

movieLength = size(s.xyzrcoPeak,3);
if isempty(whichFrames), whichFrames = 1:movieLength; end

climits = [min(s.xyzrcoPeak(4,:)) max(s.xyzrcoPeak(4,:))];

im = []; map = [];

if strcmp('birdsEye', viewPoint)
    viewPoint = [-179.6000 90];
    fnameToUse = 'ballStickBirdsEye.gif';
else
    viewPoint = [154 28];
    fnameToUse = 'ballStickSideView.gif';
end
fpathToUse = fullfile(theRefLocDir, fnameToUse);

stackInf = imfinfo(fullfile(fullfile(jlgDataDir,s.stackPath)));
stackDim.width = stackInf(1).Width;
stackDim.height = stackInf(1).Height;
stackDim.depth = length(stackInf);

saveName = [];
for ff = whichFrames
    hArr(ff) = figure('visible','off');

    hArr(ff) = xyzrcBallStickFig(s.xyzrcoPeak,ff,hArr(ff),stackDim);
    set(findobj(hArr(ff),'type','axes'),'CameraPosition', [200 200 -300.2114],...
        'View',viewPoint, 'clim', climits)
    if doShow, set(hArr(ff),'visible','on'); end
    if ff == whichFrames(end) && doSave, saveName = fpathToUse; end
    [im, map] = createGif(hArr(ff), ff, movieLength, im, map, saveName);
    pause(1)
    
    
end