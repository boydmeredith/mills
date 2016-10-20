function ml = listSubjMovieDirDates()

subj = {'J114','J115','J116','J117','J118','J122','J123'};

for ss = 1:length(subj)
    thisSubj = subj{ss};
    fldrs = dir(fullfile(jlgDataDir, thisSubj, '2015-*'));
    fldrs = fldrs([fldrs.isdir]);
    ml.(thisSubj).dates = {fldrs.name}';
    
    nDates = length(ml.(thisSubj).dates);
    
    ml.(thisSubj).hasL01Summary = false(nDates,1);
    ml.(thisSubj).hasL02Summary = false(nDates,1);
    ml.(thisSubj).hasL01Movie = false(nDates,1);
    ml.(thisSubj).hasL02Movie = false(nDates,1);
    
    for dd = 1:length(ml.(thisSubj).dates)
        thisDir = fullfile(jlgDataDir, thisSubj, ml.(thisSubj).dates{dd});
        
        ml.(thisSubj).hasL01Movie(dd) = exist(fullfile(thisDir,'L01'),'file');
        ml.(thisSubj).hasL02Movie(dd) = exist(fullfile(thisDir,'L02'),'file');
        ml.(thisSubj).hasL01Summary(dd) = exist(fullfile(thisDir,...
            'L01_referenceLocalization','summary.mat'),'file');
        ml.(thisSubj).hasL02Movie(dd) = exist(fullfile(thisDir,...
            'L02_referenceLocalization','summary.mat'),'file');
        
    end
end

