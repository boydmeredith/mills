edit function [newH, newW] = dimAfterRotation(origH, origW, angle)
% [height width] = dimAfterRotation(origH, origW, angle)
% return the new height and width of a rectangle after rotating it by the
% specified angle
%
    blockDiag = sqrt(origW^2 + origH^2);
    origAng   = asind(origH/blockDiag);
    newH      = sind(origAng + angle) * blockDiag;
    newW      = cosd(origAng - angle) * blockDiag;