function h = ballStickGifWrapper(subj, movieDate, location, whichFrames, viewNames, doSave, doShow)
% function hArr = ballStickGifWrapper(subj, movieDate, location, whichFrames, viewNames, doSave, doShow)
% 
% Plot blockwise alignments for the specified frames from a registered
% moviewith each blocks' center shown as a point in x,y,z space, with color
% corresponding to rotation and correlation magnitude shown as point size.
%
% viewNames - optional argument. Shows bird's eye view if 'birdsEye', shows
% side view otherwise.
%
%
% Saves a gif if desired



if isempty(location), location = 'L01'; end
if isempty(doSave), doSave = false; end

theRefLocDir = referenceLocalizationDir(subj, movieDate, location);
s = load(fullfile(theRefLocDir, 'summary.mat'),'xyzrcoPeak','params','stackPath');

movieLength = size(s.xyzrcoPeak,3);
if isempty(whichFrames), whichFrames = 1:movieLength; end

climits = [min(s.xyzrcoPeak(4,:)) max(s.xyzrcoPeak(4,:))];

map = [];
viewPoints = struct('birdsEye',[-179.6000 90],'sideView',[156.7 25.6]);

stackInf = imfinfo(fullfile(fullfile(jlgDataDir,s.stackPath)));
stackDim.width = stackInf(1).Width;
stackDim.height = stackInf(1).Height;
stackDim.depth = length(stackInf);

h = figure('visible','off','paperpositionmode','auto','position',[ 0 0 1000 1000] );

for ff = whichFrames
    
    h = xyzrcBallStickFig(s.xyzrcoPeak,ff,h,stackDim,s.params.mByNBlocks);
    
    for vv = 1:length(viewNames)
        thisViewName = viewNames{vv};
        set(findobj(h,'type','axes'),'View',viewPoints.(thisViewName))
        if doShow, set(h,'visible','on'); end
        [~, map] = createGif(h, ff, map, fullfile(theRefLocDir, ['ballStick_' thisViewName '.gif']));
        pause(.5)
        
    end
    delete(findall(findall(h,'Type','axe'),'Type','text'))

end