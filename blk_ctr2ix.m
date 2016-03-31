function [blk_yix blk_xix] = blk_ctr2ix(blk_cntrs, bi, blksz, imgx, imgy);


    blk_yix = [(blk_cntrs.y(bi) - blksz/2):(blk_cntrs.y(bi) + blksz/2)];
    blk_xix = [(blk_cntrs.x(bi) - blksz/2):(blk_cntrs.x(bi) + blksz/2)];

    % crop the right and bottom blocks if necessary
    % to avoiding trying to access out of bounds indices
    blk_yix = blk_yix(blk_yix <= imgy);
    blk_xix = blk_xix(blk_xix <= imgx);
    
    
end