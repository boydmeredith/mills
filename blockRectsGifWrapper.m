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

close all;

if isempty(location), location = 'L01'; end
if isempty(doSave), doSave = false; end
if ~doShow || isempty(doShow) , doShow = 'off'; else doShow = 'on'; end

theRefLocDir = referenceLocalizationDir(subj, movieDate, location);
s = load(fullfile(theRefLocDir, 'summary.mat'),'xyzrcoPeak','params','stackPath','blockLocations');
movieLength = size(s.xyzrcoPeak,3);
fullStackPath = fullfile(jlgDataDir,s.stackPath);
stackInf    = imfinfo(fullStackPath);
if isempty(whichFrames), whichFrames = 1:movieLength; end

midStackSlice = imread(fullStackPath,floor(length(stackInf)/2));

h = figure('visible',doShow,'position',[0 0 900 1000]);

map = [];

if doSave, saveName = fullfile(theRefLocDir, 'rotatingRects.gif');
else, saveName = '';
end

nBlocks = size(s.xyzrcoPeak,2);

thisCmap = colormapRedBlue;
thisCmap = thisCmap(linspace(1,size(colormapRedBlue,1),length(stackInf)),:);

bgax = axes('parent',h);
image((normalizeToZeroOne(double(midStackSlice),[0 99.5],true).^1.5)*64,'parent',bgax);
%image((normalizeToZeroOne(double(midStackSlice),5,true).^.5)*64,'parent',bgax)%,quantile(double(midStackSlice(:)),[0.0001 1-.0001]))
xlabel(bgax,'x position');ylabel(bgax,'y position');
colormap(bgax,bone);
bgaxPos = get(bgax,'position');

ax = axes('position',get(bgax,'position'),'color','none');
linkaxes([bgax ax]);

cbax = axes('parent',h,'position',[.875 bgaxPos(2) .875 bgaxPos(4)],'color','none');
axis(cbax,'off');
colormap(cbax,colormapRedBlue);
c = colorbar(cbax,'west'); caxis(cbax,[1 length(stackInf)]);

if doSave, fprintf('\n\tWriting images...');end
tic
for thisFrameNo=whichFrames,
    %%

    
    title(ax,sprintf('%s %s frame: %03d/%03d',subj,movieDate,thisFrameNo,whichFrames(end)));
    for thisBlockNo = 1:nBlocks,
        thisColor=thisCmap(s.xyzrcoPeak(3,thisBlockNo,thisFrameNo),:);
        bInf=getBlockInf(s.blockLocations(:,:,thisBlockNo));
        
        drawShadedRect(bInf.width,bInf.height,s.xyzrcoPeak(1,thisBlockNo,thisFrameNo),...
            s.xyzrcoPeak(2,thisBlockNo,thisFrameNo),s.xyzrcoPeak(4,thisBlockNo,thisFrameNo),...
            thisColor,s.xyzrcoPeak(6,thisBlockNo,thisFrameNo), ax);
        %drawnow

    end;
    axis([bgax ax],'image');
    padWidth = 50;
    set( ax,'ylim',[1 size(midStackSlice,1)]+[-padWidth padWidth],'xlim',[1 size(midStackSlice,2)]+[-padWidth padWidth])

    %%
    pause(.55);
    set(c,'ydir','rev')
    [~, map] = createGif(h, thisFrameNo, map, saveName);

%%
    delete(findall(findall(ax,'Type','axe'),'Type','text'))
    cla(ax);
    
    if ~mod(thisFrameNo,15)
        toc
        pause(1)
        tic
        fprintf('%d... ',thisFrameNo);
    end
end
close all; 


