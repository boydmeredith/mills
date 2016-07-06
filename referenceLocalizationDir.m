function thePath = referenceLocalizationDir(subj, date, location)
if nargin < 3, location = 'L01'; end
thePath = fullfile(jlgDataDir,subj,date,sprintf('%s_referenceLocalization',location));
end