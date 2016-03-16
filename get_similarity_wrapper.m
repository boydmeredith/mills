function [] = get_similarity_wrapper(subj)
% get_similarity_wrapper(subj)
%
% 
%
% input
% subj       subject ID. eg 'J115' 
%

setup; % add useful paths; configure default figure text

off_range = -50:50; % range of offsets to search for block matches
blksz     = 100;    % size of blocks to divide images into
z_range   = -3:3;   % range around best matches to use

% load the z stack into a 3d array
stack_names = dir(fullfile(data_dir,subj,'2015-*-*__post_stack_00*_AVERAGE.tif'));
stack_date  = stack_names(1).name(1:10);
stack_path  = fullfile(data_dir, subj, stack_names(1).name);
stack       = load_stack(stack_path);

% divide into blocks (define centers in reference to the stack images)
stack_info  = imfinfo(stack_path);
blk_cntrs   = make_blocks(stack_info(1).Height, stack_info(1).Width, blksz);
nblks       = length(blk_cntrs.x);

% get a list of movies each corresponding to a day of imaging
movie_names = dir(fullfile(data_dir,subj,['2015-*__L01__AVERAGE.tif']));
movie_names = {movie_names(:).name};
ndays       = length(movie_names);

% pre-allocate best z choice array
blks_best_z = nan(nblks,ndays);

for di = 1:ndays
    % pre-allocate correlations array 
    C = nan(nblks, length(off_range), length(off_range), length(stack_info));
    
    % load first frame of this day's movie
    img_path = fullfile(data_dir,subj,movie_names{di});
    img_name = movie_names{di};
    img_date = img_name(1:10);
    img = load_and_norm_img(img_path);
    [imgy imgx] = size(img);

    % let us know what we're image we're working on
    fprintf('\n\nworking on %s...', img_name);

    % if it's the first day, find the best match in the stack and 
    % guess that it's the best z for all blocks
    if di == 1
        best_z = match_to_stack(stack, img);
        blks_best_z(:,1) = best_z;
    end

    % get cross correlation for each block
    for bi = 1:nblks
        fprintf('\n\tblock %i...', bi);
        % create block indices around this block center
        blk_yix = [(blk_cntrs.y(bi) - blksz/2):(blk_cntrs.y(bi) + blksz/2)];
        blk_xix = [(blk_cntrs.x(bi) - blksz/2):(blk_cntrs.x(bi) + blksz/2)];
       
        % crop the right and bottom blocks if necessary
        % to avoiding trying to access out of bounds indices
        blk_yix = blk_yix(blk_yix <= imgy);
        blk_xix = blk_xix(blk_xix <= imgx);

        % select block
        block   = img(blk_yix,blk_xix);
        
        % determine the neighborhood of z-values to look at
        z_to_check = blks_best_z(bi,di) + z_range;
        % get the cross correlation of this block against the stack image
        [c xoff yoff]   = get_xyz_similarity(stack(:,:,z_to_check), block, blk_cntrs.x(bi), blk_cntrs.y(bi), off_range);

        C(bi,:,:,z_to_check) = c;

        % update best z with the z that has the max correlation
        [~, maxix] = max(c(:));
        [~, ~, zix] = ind2sub(size(c),maxix);
        fprintf('\n\tfavorite z for this block in neighborhood around %i...', blks_best_z(bi, di));
        blks_best_z(bi,di) = z_to_check(zix); 
        if bi < nblks
            blks_best_z(bi,di+1) = blks_best_z(bi,di);
        end
        fprintf('s %i', blks_best_z(bi, di));

    end
    % now let's save this correlation matrix
    savepath = fullfile(data_dir,subj,[img_name(1:end-4) '_correlations.mat']); 
    save(savepath, 'C','off_range')
end


end
