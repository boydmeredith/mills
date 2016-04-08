function [C yoff xoff] = get_xyz_similarity(stack,blk_img, blk_cntr_x, blk_cntr_y, off_range)
% [C xoff yoff] = get_xyz_similarity(stack,blk_img, blk_cntr_x, blk_cntr_y, off_range)
%  
% find the cross correlation within some range of allowable offsets 
% between a small block image and slices from a zstack. 
% 
% input
% stack         the slices that the block should be compared to
% blk_img       the small block from the image we are trying to register
% blk_cntr_x    x coordinate in stack to calculate (and restrict) offsets 
% blk_cntr_y    y coordinate in stack to calculate (and restrict) offsets 
% off_range     allowed x,y offsets for maximizing cross correlation


n_z = size(stack,3);
% pre-allocate correlations array 
C = nan(length(off_range), length(off_range), n_z);

for zx = 1:n_z
    
        stack_slice = stack(:,:,zx);
        
        % get the cross correlation of this block against the stack image
        [c yoff xoff]   = block_xy_similarity(stack_slice, blk_img, blk_cntr_x, blk_cntr_y);

        % store xcorr values for allowable offsets
        xix = find(ismember(xoff,off_range));
        yix = find(ismember(yoff,off_range));

        % pad c with nans so that it fits as expected even if we've cropped the blocks
        c = padarray(c(yix,xix),[length(off_range)-length(yix) length(off_range)-length(xix)],nan,'post');

        C(:,:,zx) = c;
end

end
