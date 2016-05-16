function [yy, xx]  = getBlockIx(blockLocBin)
% [yy xx] = getBlockIx(blockLocBin)
% find max across all columns and all rows of a binary matrix to get
% indices of the nonzero entries
    yy = find(max(blockLocBin,[],2));
    xx = find(max(blockLocBin,[],1));