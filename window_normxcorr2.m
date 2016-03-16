function [maxc d]  = window_normxcorr2(temp,slice_img,max_d)
    
    if numel(temp) < numel(slice_img)
        c = (normxcorr2(temp,slice_img));
        smallx = size(temp,2); smally = size(temp,1);
    elseif numel(temp) > numel(slice_img)
        c = (normxcorr2(slice_img,temp));
        smallx = size(slice_img,2); smally = size(slice_img,1);
    else
        c = (normxcorr2(temp([1 size(temp,1)-1],[1 size(temp,2)-1]),...
            slice_img));
        smallx = size(temp,2)-1; smally = size(temp,1)-1;
    end
    
    yoffsets = [1:size(c,1)]- smally;
    xoffsets = [1:size(c,2)]- smallx;
    [X Y ] = meshgrid(xoffsets, yoffsets);
    c(X.^2 + Y.^2 > max_d^2) = -Inf;
    
    maxc              = max(c(:));
    [d.ypeak d.xpeak] = find(c==maxc);
    d.d = Inf;
    for i = 1:length(d.ypeak)
        
        this_yoff            = d.ypeak(i) - smally;
        this_xoff            = d.xpeak(i) - smallx;
        this_d               = sqrt(this_yoff^2 + this_xoff^2);
        if this_d < d.d
            d.yoff = this_yoff;
            d.xoff = this_xoff;
            d.d    = this_d;
        end
    
end