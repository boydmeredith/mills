function [b]  = getBlockInf(blockLoc)
% [b] = getBlockIInf(blockLocBin) 
% Takes a logical matrix specifying a block's location and returns a struct
% containing information about a blocks location with fields:
%   indY        - indices in Y dimension
%   indX        - it's indices in the X dimension
%   height      - size in Y dimension
%   width       - size in X dimension
%   ctrY        - Y coordinate of block center
%   ctrX        - X coordinate of block center
%   blockLocBin - the logical matrix provided as input to the function
    b.indX     = find(max(blockLoc,[],1));
    b.indY     = find(max(blockLoc,[],2));
    b.width    = length(b.indX);
    b.height   = length(b.indY);
    b.ctrX     = mean(b.indX);
    b.ctrY     = mean(b.indY);
    b.blockLoc = blockLoc;