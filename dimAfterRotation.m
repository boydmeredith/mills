function [height width] = dimAfterRotation(blockH, blockW, angle)
    blockDiag = sqrt(blockW^2 + blockH^2);
    origAng   = asind(blockH/blockDiag);
    height    = sind(origAng + angle) * blockDiag;
    width     = cosd(origAng - angle) * blockDiag;