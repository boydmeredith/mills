function plotAllBlockCenters(thisSubj)

dateList = listSubjMovieDirDates
dateList = dateList.(thisSubj).dates;
f = figure;
h = makeSubplots(f,5,ceil(length(dateList)/5),.01, .01, [0 0 .95 .95]);




for dd = 1:length(dateList)
    thisDate = dateList{dd};
     
    location = '';
    
    if exist(fullfile(referenceLocalizationDir(thisSubj,thisDate,'L01'),...
            'summary.mat'), 'file');
        location = 'L01';
        plotBlockCenters(thisSubj, thisDate, location, 'ax', h(dd));
    elseif exist(fullfile(referenceLocalizationDir(thisSubj,thisDate,'L02'),...
            'summary.mat'), 'file');
        location = 'L02';
        plotBlockCenters(thisSubj, thisDate, location, 'ax', h(dd));
    end

    
    title(h(dd), [thisSubj ' - ' thisDate ' - ' location]);
end

