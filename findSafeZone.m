%function safeImg = findSafeZone(shifts, meanImg)
   % safeImg = findSafeZone(shifts, meanImg)
   % trim max shift from the corners of the meanImg. Used to return
   % the mean of a movie within a zone that is always guaranteed to be
   % in the frame (according to motion correction)
    
   
%    
%    maxShifts = ceil(max(abs(shifts)));   
%    [h w] = size(meanImg);
%    safeImg = meanImg([maxShifts(1):(h-maxShifts(1))],...
%             [maxShifts(2):(w-maxShifts(2))]);




%%
uncorrectedMoviePath = '/Volumes/tank/jlgauthi/Data/J115/2015-09-25/L01__UNCORRECTED_AVERAGE.gif';
mcDir = '/Volumes/tank/jlgauthi/Data/J115/2015-09-25/L01/';

uncorrectedMovie = squeeze(imread(uncorrectedMoviePath, 'Frames', 'all'));

nFrameAv = 1000;
nAvFramesInMovie  = size(uncorrectedMovie,3); 
allShifts = zeros(nAvFramesInMovie,nFrameAv,2);
%frameShifts = nan(length(uncorrectedMovie),2);
mcFiles = dir(fullfile(mcDir,'ac*mat'));
mcFiles = {mcFiles.name};

%%     
% use currentFrame to keep track of where we are in motion corrected movie
% which has 1000x as many frames as the averaged movie (or in some cases
% 3000x)
[h w t] = size(uncorrectedMovie);
safeZone = false(h,w,t);
currentFrame = 1;
for ff = 1:numel(mcFiles)
    mc = load(fullfile(mcDir,mcFiles{ff}), 'shifts');
    nFrames = size(mc.shifts,1);
    nAvFrames = ceil(nFrames/nFrameAv);
    for subf = 1:nAvFrames
        first = nFrameAv*(subf-1)+1;
        last  = min(nFrames, first+nFrameAv-1);
        thisNFrames = last - first + 1;
        allShifts(currentFrame, 1:thisNFrames, :) = ...
            mc.shifts(first:last,:);
        
        maxYShift = max( ceil(max(allShifts(currentFrame,:,1),[],2)), 1);
        minYShift = min( floor(min(allShifts(currentFrame,:,1),[],2)), h);

        maxXShift = max( ceil(max(allShifts(currentFrame,:,2),[],2)), 1);
        minXShift = min( floor(min(allShifts(currentFrame,:,2),[],2)), w);
        
        safeZone([maxYShift:(h+minYShift)],[maxXShift:(w+minXShift)], currentFrame) = true;
        
        currentFrame = currentFrame+1;

    end
end









