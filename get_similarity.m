function [C theta_scores] = get_similarity(stack_path,img_path, z_range)
% get_similarity(stack_path,img_path,z_range)
% 
%
%
%
%
%
% allowed x,y offsets for maximizing cross correlation
xoff_range  = -50:50;
yoff_range  = -50:50;
theta_range = -20:20;
blksz       = 100;

% TODO: need to choose z range, based on previous information

% load image
img = load_and_norm_img(img_path);
% get image dimensions
[imgy imgx] = size(img);

% divide into blocks (define centers in reference to the stack images)
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

        % create block indices 
        blk_yix = [(blk_cntrs.y(bx) - blksz/2):(blk_cntrs.y(bx) + blksz/2)];
        blk_xix = [(blk_cntrs.x(bx) - blksz/2):(blk_cntrs.x(bx) + blksz/2)];
       
        % crop the right and bottom blocks to avoiding trying to access out of bounds indices
        blk_yix = blk_yix(blk_yix <= imgy);
        blk_xix = blk_xix(blk_xix <= imgx);

        % select block
        block   = img(blk_yix,blk_xix);

        figure(1); clf; imshow(block);  
        % get the cross correlation of this block against the stack image
        [c xoff yoff]   = block_xy_similarity(stack_slice, block, blk_cntrs.x(bx), blk_cntrs.y(bx));

        % store xcorr values for allowable offsets
        xix = find(ismember(xoff,xoff_range));
        yix = find(ismember(yoff,yoff_range));

        % pad c with nans so that it fits as expected even if we've cropped the blocks
        c = padarray(c(yix,xix),[length(yoff_range)-length(yix) length(xoff_range)-length(xix)],nan,'post');

        C(bx,:,:,zx) = c;
        keyboard;
end

end
