if exist('/Volumes/jlgauthi','dir')
    addpath(genpath('/Volumes/jlgauthi/code'));
elseif exist('/usr/people/jlgauthi/code','dir')
    addpath(genpath('/usr/people/jlgauthi/code'));
else
    fprintf('ERROR: directory not mounted\n');
end
