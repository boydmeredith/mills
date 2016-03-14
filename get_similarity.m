function [C theta_scores] = get_similarity(stack_path,img_path, z_range)
% allowed x,y offsets for maximizing cross correlation
xoff_range  = -50:50;
yoff_range  = -50:50;
theta_range = -20:20;
blksz       = 100;

% TODO: need to choose z range, based on previous information

% load image
img = load_and_norm_img(img_path);

% divide image into blocks (define centers in reference to the stack images)
stack_info  = imfinfo(stack_path);
blk_cntrs   = make_blocks(stack_info(1).Height, stack_info(1).Width, blksz);
nblks       = length(blk_cntrs.x);

% pre-allocate correlations array 
C = nan(nblks, length(xoff_range), length(yoff_range), length(stack_info));

for zx = z_range
    stack_slice     = load_and_norm_img(stack_path,zx);
    for bx = 1:nblks
        % TODO: should check to see if we need to check zstack slice zx for this block
        % TODO: add rotation

        % get the indices of the current block
        blk_yix = [(blk_cntrs.y(bx) - blksz/2):(blk_cntrs.y(bx) + blksz/2)];
        blk_xix = [(blk_cntrs.x(bx) - blksz/2):(blk_cntrs.x(bx) + blksz/2)];
        % TODO: fix index exceeding matrix dimension 
        block   = img(blk_yix,blk_xix);
        figure(1); clf; imshow(block);  
        % get the cross correlation of this block against the stack image
        [c xoff yoff]   = get_xy_similarity(stack_slice, block);
        % correct the offsets returned 
        xoff = xoff - max(blk_xix);
        yoff = yoff - max(blk_yix);
        
        % store xcorr values for allowable offsets
        xix = find(ismember(xoff,xoff_range));
        yix = find(ismember(yoff,yoff_range));
        C(bx,:,:,zx) = c(yix,xix);
        pause(.5);
end

end
