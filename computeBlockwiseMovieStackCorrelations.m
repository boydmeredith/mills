function [] = computeBlockwiseMovieStackCorrs(moviePath, stackPath, varargin)

    corrType = 'uint8';

    % things that should be in varargin
    nBlockSpan = 10;
    overlap    = .2;
    rotationAngles  = -10:.25:10;
    rZero = find(rotationAngles == 0);
    nRotation = length(rotationAngles);
    nBlocks = nBlockSpan^2;
    zNeighborhood
    rNeighborhood
    xyNeighborhood
    
    % load stack (expect gif)
    stack = imread(stackPath,'Frames','all');
    
    % crop dark edge off stack by removing rows of all zero entries
    stack = cropStack(stack);
    
    % get stack size
    [stackHeight, stackWidth, stackDepth] = size(stack);    
    
    % load movie
    movie = imread(moviePath,'Frames','all');
    % get movie size
    [movieHeight, movieWidth, movieLength] = size(movie);
    
    % divide movie into blocks by specifying the dimensions of the movie,
    % the number of vertical and horizontal bars to use to make the blocks,
    % the percent overlap the bars should have, and the amount of rotation
    % that we should be able to accomodate at the margins of the image
    blockLocations = makeBlockLocations(movieHeight, movieWidth, ...
                                        nBlockSpan, overlap, max(rotationAngles));
    
    % initialize matrices to hold correlations and the xyzr values that
    % maximize those correlations
    corrMat      = zeros(movieHeight,movieWidth,stackHeight,nRotation,...
                        movieLength, nBlocks, corrType);
                    
    xyzrLocation = zeros(4,movieLength,nBlocks);
    
    % iterate through the frames of the movie 
    
    for ff = 1:movieLength,
        % select the relevant movie frame
        movieFrame = movie(:,:,ff);
        
        for bb = 1:nBlocks,
            % select the relevant block location matrix
            blockLoc = blockLocations(:,:,bb);
            
            for zz = 1:stackDepth,
                % select the relevant stack slice
                stackSlice = stack(:,:,zz);
                
                % compute correlations for this block without applying
                % rotation and store in correlation matrix
                corrMat(:,:,zz,rZero,ff,bb) = computeBlockImageCorrs(movieFrame, ...
                    blockLoc, 0, stackSlice, corrType);
                
            end
            
            % find the z that best matches the unrotated block and then
            % attempt registration within that block by
            % using rotation 
            [peaksNoRot] = getPeakCorrCoords(corrMat(:,:,:,rZero,ff,bb), blksz, rotationAngles);
            zPeakNoRot = peaksNoRot(3);
            bestStackSliceNoRot = stackSlice(:,:,zPeakNoRot)
            
            for rr = 1:length(rotationAngles)
                % rotate the block by rotAngle and put into corrMat at the
                % rotation indexed by rr
                rotAngle = rotationAngles(rr);
                blockRot = rotateAndReselectBlock(movieFrame, blockLoc, rotAngle);
                corrMat(:,:,zz,rr,ff,bb) = computeBlockImageCorrs(movieFrame, ...
                    blockLoc, rotAngle, stackSlice, corrType);

            end
            
            for zz = neighborhoodAroundBestZ
                for rr = neighborhoodAroundBestR
                    % select region of stack slice in xy based on how much
                    % we want to fill in the corrMat. throw away
                    % appropriate indices and put into appropriate of
                    % persistent corrMat
                    
                end
            end
            
            xyzrLocations(:,ff,bb) = getPeakCorrCoords(corrMat(:,:,:,:,ff,bb), blksz, rotationAngles);
                
        end
    end
    
    
%     % to save
%     blockLocations
%     corrMat
%     xyzLocation (should correspond to the center of the block)
%     rotationAngles
%     stackPath
%     analysisDate

end

function [peaks yPeakInRef xPeakInRef rPeakInRef] = getPeakCorrCoords(corrMat, blksz, rotationAngles)
    peaks = ind2sub(size(corrMat),find(corrMat==max(corrMat(:))));
    yPeakInRef = yPeak - blksz + 1;
    xPeakInRef = xPeak - blksz + 1;
    rPeakInRef = rotationAngles(rPeak);
end

    