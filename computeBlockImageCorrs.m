function [corrMat] = computeBlockImageCorrs(block, reference, corrType);
% corrMat = computeBlockImageCorrs(block, reference, corrType); compute the
% normalized 2d cross correlation between a block and a reference. The
% block is expected to be smaller than the reference. 
% 
% inputs:
% block      -      An image selected from a movie frame that we want to
%                   register to the reference
% reference  -      An image in which to look for a match of the block
% corrType   -      A variable type for storage of the correlation matrix


corrMat = normxcorr2(block,reference);
if ~strcmp(corrType,'double')
    corrMat = cast(double(intmax(corrType)) * corrMat,corrType);
end

end