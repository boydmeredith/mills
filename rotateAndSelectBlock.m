function [block] = rotateAndSelectBlock(movieFrame, blockLoc, angle);
% block = rotateAndReselectBlock(movieFrame, blockLoc, angle) 
% 
% Rotate movieFrame clockwise about the center of the blockLoc by the
% specified angle and return the new image in the blockLoc
%
% movieFrame -          An image from which to select the block 
% blockLoc   -          A logical matrix of the same size as movieFrame with
%                       1's indicating the position of the bloc
% angle      -          An angle in degrees by which to rotate the image
%
% Note: assumes that margins have been built into the edges of the
% blockLocations so that the rotation will not introduce empty pixels to
% the blockT = maketform('affine',[cosd(r) sind(r) 0; -sind(r) cosd(r) 0; 0 0 1]);

% get dimensions of the movie frame
[movieFrameHeight movieFrameWidth] = size(movieFrame);

% ensure the movie frame and block location matrix have same dimensions
assert(isequal(size(blockLoc), size(movieFrame)));

% get the indices, dimensions, and center of the block in the movie Frame
bInf = getBlockInf(blockLoc);

T = maketform('affine',[cosd(angle) sind(angle) 0; -sind(angle) cosd(angle) 0; 0 0 1]);


block =imtransform(movieFrame,T,'bilinear',...
    'udata',[1 movieFrameWidth]  - bInf.ctrX,...
    'vdata',[1 movieFrameHeight] - bInf.ctrY,...
    'xdata',([1 bInf.width]) -(bInf.width+1)/2,...
    'ydata',([1 bInf.height])-(bInf.height+1)/2,...
    'xyscale',1);