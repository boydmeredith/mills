function [im, map]  = createGif(fig, frameNo, map, gifName)

% note on hardcopy: matlab says not to use it, but the alternative getframe
% requires you to have the figure open and fully visible to the user.
% the other alternative is print, but I can't figure how to collect
% its result in a variable
if ~ishandle(fig)
    fig = figure(fig,'Visible','off');
end

f.cdata = hardcopy(fig, '-Dopengl', '-r0');
if frameNo == 1
    [im, map] = rgb2ind(f.cdata,256,'nodither');
    if ~isempty(gifName)
        imwrite(im,map,gifName,'gif','DelayTime',.2,'LoopCount',Inf);
    end
else
    im(:,:,1) = rgb2ind(f.cdata,map,'nodither');
    if ~isempty(gifName)
        imwrite(im,map,gifName,'gif','DelayTime',.2,'WriteMode','append');
    end
end

end