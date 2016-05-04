function safeImg = findSafeZone(shifts, meanImg)
   % safeImg = findSafeZone(shifts, meanImg)
   % trim max shift from the corners of the meanImg. Used to return
   % the mean of a movie within a zone that is always guaranteed to be
   % in the frame (according to motion correction)
    
   maxShifts = ceil(max(abs(shifts)));   
   [h w] = size(meanImg);
   safeImg = meanImg([maxShifts(1):(h-maxShifts(1))],...
            [maxShifts(2):(w-maxShifts(2))]);





