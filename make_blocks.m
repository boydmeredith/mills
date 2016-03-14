function [pts] = make_blocks( h, w, blcksz, center)
if nargin < 4
    center =1 ;
end
% function [pts.x pts.y] = make_blocks(w, h, blcksz)
% figure out how to evenly tile a rectangle with square blocks with sides
% blcksz

n_blocks_x = floor(w / blcksz);
n_blocks_y = floor(h / blcksz);

start_x = floor(mod(w,blcksz)/2);
start_y = floor(mod(h,blcksz)/2);

if center
    [pts.x pts.y] = meshgrid([start_x+blcksz/2:blcksz:w],[start_y+blcksz/2:blcksz:h]);
else
    [pts.x pts.y] = meshgrid([start_x:blcksz:w-blcksz],[start_y:blcksz:h-blcksz]);
end
pts.x = reshape(pts.x,1,numel(pts.x));
pts.y = reshape(pts.y,1,numel(pts.y));



%scatter(pts.x,pts.y,'k.')
