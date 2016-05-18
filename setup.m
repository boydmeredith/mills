addpath(genpath('~/tank/mills/'));
jlgauthi_util_dir = '/Volumes/tank/jlgauthi/code/utility';
if exist(jlgauthi_util_dir)
    addpath(jlgauthi_util_dir);
else
    addpath('utility');
end
data_dir = '~/Documents/Princeton/Rotations/tank/';
data_dir = '/Volumes/tank/jtb3/Data/';
easy_subject = 'J115';
hard_subject = 'J114';

set(0,'defaultaxesfontsize',15);
set(0,'defaulttextfontname','helvetica');
