function [C] = get_xy_similarity(ref,img)
% C = get_xy_similarity(ref,img)
% 
% get the 2d cross correlation between two images as well the 
% x and y offsets the correspond to each cell in the xcorr matrix 
%
% input:
%
% ref       a slice from the stack to compare the img to
% img       a slice from a day of experiment to match the reference
%
% output:
% C         a matrix of correlations for a variety of x and y offsets
% note that max cross correlation corresponds to the estimated 
% location of the lower-right corner of the img

% mean center both images
ref = ref - mean(ref(:));
img = img - mean(img(:));

% get the cross correlation between the images
C = xcorr2(ref,img);

end
