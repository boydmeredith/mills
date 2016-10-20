function edgeBlocks = getEdgeBlocks(gridSize)
edgeBlocks = logical(reshape(padarray(zeros(gridSize-[2 2]),[1 1],1,'both'),1,prod(gridSize)));
end