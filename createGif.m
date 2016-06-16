function [im, map]  = createGif(fig, frameNo, movieLength, im, map, gifName)

% note on hardcopy: matlab says not to use it, but the alternative getframe
% requires you to have the figure open and fully visible to the user.
% the other alternative is print, but I can't figure how to collect
% its result in a variable
if ~ishandle(fig)
    fig = figure(fig,'Visible','off');
end

f.cdata = hardcopy(fig, '-Dzbuffer', '-r0');
if frameNo == 1
    [im, map] = rgb2ind(f.cdata,256,'nodither');
    im(1,1,1,movieLength) = 0;
else
    im(:,:,1,frameNo) = rgb2ind(f.cdata,map,'nodither');
end
if ~isempty(gifName)
    imwrite(im,map,gifName,'DelayTime',.2,'LoopCount',Inf);
end
end