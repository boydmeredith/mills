function [ii,jj] = findClosestAcceptablePoint(initialPoint,acceptableZone)

% function [ii,jj] = findClosestAcceptablePoint(initialPoint,acceptableZone)
% Given an initial point (i,j), find the closest acceptable point in a space as
% specified by the logical matrix acceptableZone 
%

[nY,nX] = size(acceptableZone);

if acceptableZone(round(initialPoint(1)),round(initialPoint(2)))
    ii = initialPoint(1);
    jj = initialPoint(2);
    return
end


% get distance to every point
% this line is slow:
% dists = sqrt(bsxfun(@plus,((1:nX)-initialPoint(2)).^2,((1:nY)-initialPoint(1))'.^2));
% this line is half as slow:
dists = sqrt(repmat(((1:nX)-initialPoint(2)).^2,nY,1)+repmat(((1:nY)-initialPoint(1))'.^2,1,nX));

% set unacceptable points to have distance of inf
dists(~acceptableZone) = inf;


[~,closestPoint] = min(dists(:));

[ii,jj] = ind2sub([nY nX],closestPoint);

%[ii,jj] = find(dists == min(min(dists)));


