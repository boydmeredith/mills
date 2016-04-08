function [C yoff xoff] = block_xy_similarity(ref, blk_img, cntr_x, cntr_y)
%  [C yoff xoff] = block_xy_similarity(ref, blk_img, blk_cntr)
% 
% get the 2d cross correlation of a small image and a larger reference
% image. return offsets relative to an x, y position in the reference
%
% input:
% ref       large reference image (maybe a slice from a zstack?)
% blk_img   small image block to situate within the ref image
% cntr_x     x position to calculate offsets with respect to
% cntr_y     y position to calcualte offsets with respect to
%
% output:
% C         matrix of cross correlations
% yoff      array of offsets from the blk_cntrs y coordinate
% xoff      array of offsets from the blk_cntrs x coordinate

% get the dimensions of each image
[refy refx] = size(ref);
[blky blkx] = size(blk_img);

% get bottom right corner of block
blk_br_x = cntr_x + floor(blkx/2);
blk_br_y = cntr_y + floor(blky/2);
C = normxcorr2(blk_img, ref);

% figure out what offsets these correspond to
% the indices of C correspond to positions of the bottom right corner
% relative to the reference (so, at 1,1 the block overlaps with the
% reference only by its bottom left pixel)
[Cy Cx] = size(C); 

% subtract the bottom right corner of the reference block center in
% order to get the offsets in terms of that that reference block
% position
yoff = [1:Cy]-blk_br_y;
xoff = [1:Cx]-blk_br_x;


end
