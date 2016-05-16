function stack = cropStack(stack)
% stack = cropStack(stack)
% looks for rows or columns in the image that
% contain nothing and gets rid of them


stackMax = max(squeeze(stack),[],3);

[y x] = size(stackMax);

minNanRow = y+1;
minNanCol = x+1;

emptyRows = ismember(stackMax, zeros(1,x),'rows');
emptyCols = ismember(stackMax', zeros(1,y),'rows');
if min(find(emptyRows)) < minNanRow
    minNanRow = min(find(emptyRows));
end
if min(find(emptyCols)) < minNanCol
    minNanCol = min(find(emptyCols));
end
       
stack = stack(1:(minNanRow-1),1:(minNanCol-1),:);

end
