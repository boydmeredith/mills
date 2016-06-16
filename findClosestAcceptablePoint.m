function [ii,jj] = findClosestAcceptablePoint(initialPoint,acceptableZone)

% function [ii,jj] = findClosestAcceptablePoint(initialPoint,acceptableZone)
% Given an initial point (i,j), find the closest acceptable point in a space as
% specified by the logical matrix acceptableZone 
%

[nY,nX] = size(acceptableZone);

% get distance to every point
dists = sqrt(bsxfun(@plus,((1:nX)-initialPoint(2)).^2,((1:nY)-initialPoint(1))'.^2));

% set unacceptable points to have distance of inf
dists(~acceptableZone) = inf;


[~,closestPoint] = min(dists(:));

[ii,jj] = ind2sub([nY nX],closestPoint);

%[ii,jj] = find(dists == min(min(dists)));


