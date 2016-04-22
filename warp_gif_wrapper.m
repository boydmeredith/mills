function [im] = warp_gif_wrapper(subj, day)

setup;
if ~isempty(day)
    cmat_dir = fullfile(data_dir, subj, ['2015-' day '__L01__AVERAGE_corrs']);
    cmat_names = dir(fullfile(cmat_dir, '*frame*mat'));
    cmat_names = {cmat_names(:).name};
    gifname = 'warping_within_day';
else
    cmat_dir = fullfile(data_dir, subj);
    cmat_day_dirs = dir(fullfile(cmat_dir, ['2015-*__L01__AVERAGE_corrs']));
    cmat_day_dirs(~[cmat_day_dirs.isdir]) = [];
    cmat_names = {};
    for k = 1:length(cell(length(cmat_day_dirs)))
        thisdir = cmat_day_dirs(k).name;
        f = dir(fullfile(cmat_dir,thisdir,'*frame001.mat')); 
        cmat_names = [cmat_names strcat([thisdir filesep],{f.name})];
    end
    gifname = 'warping_x_days'
end
make_warp_gif(cmat_dir, cmat_names, gifname);
