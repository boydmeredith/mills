function [block] = rotateAndSelectBlock(movieFrame, bInf, angle)
% block = rotateAndReselectBlock(movieFrame, bInf, angle) 
% 
% Rotate movieFrame clockwise about the center of the bInf.blockLoc by the
% specified angle and return the new image in the bInf.blockLoc
%
% movieFrame -          An image from which to select the block 
% bInf       -          Structure containing a logical matrix blockLoc of the same
%                       size as movieFrame with 1's indicating the position of the block 
%                       The struct also has the center coordinates and
%                       dimensions of the block. If user just passes in the
%                       logical blockLoc matrix, the struct will be
%                       generated.
% angle      -          An angle in degrees by which to rotate the image
%
% Note: assumes that margins have been built into the edges of the
% bInf.blockLocations so that the rotation will not introduce empty pixels to
% the blockT = maketform('affine',[cosd(r) sind(r) 0; -sind(r) cosd(r) 0; 0 0 1]);

% get dimensions of the movie frame
[movieFrameHeight movieFrameWidth] = size(movieFrame);

% get the indices, dimensions, and center of the block in the movie Frame
if islogical(bInf)
    bInf = getBlockInf(bInf);
end

% ensure the movie frame and block location matrix have same dimensions
assert(isequal(size(bInf.blockLoc), size(movieFrame)));


T = maketform('affine',[cosd(angle) sind(angle) 0; -sind(angle) cosd(angle) 0; 0 0 1]);


block =imtransform(movieFrame,T,'bilinear',...
    'udata',[1 movieFrameWidth]  - bInf.ctrX,...
    'vdata',[1 movieFrameHeight] - bInf.ctrY,...
    'xdata',([1 bInf.width]) -(bInf.width+1)/2,...
    'ydata',([1 bInf.height])-(bInf.height+1)/2,...
    'xyscale',1);