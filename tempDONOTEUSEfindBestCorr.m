function thisCorr = tempDONOTEUSEfindBestCorr(rr,movieFrame,stackSlice,thisBlockLoc,nbrhdInf,pRes)


rr = max(rr,-10);
rr = min(rr,10);



% rotate the block according to rr
blockRot = rotateAndSelectBlock(movieFrame, thisBlockLoc, rr);

% use find to get y and x indices as well as the values
% themselves
[~,~, thisCorr] = find(computeBlockImageCorrs(blockRot, ...
    stackSlice, nbrhdInf, pRes.minCorrOverlap, pRes.corrType));
thisCorr = double(thisCorr)/double(intmax('uint16'));
thisCorr = double(-max(thisCorr(:)));
disp(thisCorr)

