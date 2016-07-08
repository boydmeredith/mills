function [] = blockRectsGifWrapper(subj,movieDate,location,whichFrames,doSave,doShow)
% function [] = blockRectsGifWrapper(subj,movieDate,location,whichFrames,doSave,doShow)
% example usage blockRectsGifWrapper('J114','2015-09-18',[],true,false)
%
% Creates a gif showing the position, rotation, and z depth of each
% registered block for a given subject, movie, and recording location
% 
% doSave determines whether the movie will be saved (in its
% referenceLocalization directory)
%
% doShow determines whether the plots will be displayed as they are created


if isempty(location), location = 'L01'; end
if isempty(doSave), doSave = false; end
if ~doShow || isempty(doShow) , doShow == 'off'; else doShow = 'on'; end

theRefLocDir = referenceLocalizationDir(subj, movieDate, location);
s = load(fullfile(theRefLocDir, 'summary.mat'),'xyzrcoPeak','params','stackPath','blockLocations');
movieLength = size(s.xyzrcoPeak,3);
if isempty(whichFrames), whichFrames = 1:movieLength; end

h = figure('visible',doShow,'paperpositionmode','auto','position',[ 0 0 1000 1000]);
ax = axes('parent',h);
im = []; map = [];

fpathToUse = fullfile(theRefLocDir, 'rotatingRects.gif');

nBlocks = size(s.xyzrcoPeak,2);

for thisFrameNo=whichFrames,
    cla(ax);
    colormap(ax,colormapRedBlue);
    colorbar(ax);
    caxis(ax,[s.params.whichSlices(1) s.params.whichSlices(end)]);
    set(ax, 'ylim', [1 512], 'xlim', [1 512]);
    xlabel(ax,'x position');ylabel(ax,'y position');
    title(ax,sprintf('%s %s frame: %03d/%03d',subj,movieDate,thisFrameNo,whichFrames(end)));
    for thisBlockNo = 1:nBlocks,
        bInf=getBlockInf(s.blockLocations(:,:,thisBlockNo));
        drawShadedRect(bInf.width,bInf.height,s.xyzrcoPeak(1,thisBlockNo,thisFrameNo),...
            s.xyzrcoPeak(2,thisBlockNo,thisFrameNo),s.xyzrcoPeak(4,thisBlockNo,thisFrameNo),...
            s.xyzrcoPeak(3,thisBlockNo,thisFrameNo),s.xyzrcoPeak(6,thisBlockNo,thisFrameNo), ax);
        
    end;
    alpha(findobj(ax),.5);
    pause(.2);
    if ~mod(thisFrameNo,15) || thisFrameNo==whichFrames(end) && doSave, 
        saveName = fpathToUse; 
    else
        saveName = [];
    end
    [im, map] = createGif(h, thisFrameNo, movieLength, im, map, saveName);
end