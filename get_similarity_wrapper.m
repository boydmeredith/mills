function [] = get_similarity_wrapper(subj, days, timepoints)
% get_similarity_wrapper(subj)
%
% 
%
% input
% subj       subject ID. eg 'J115' 
%
setup; % add useful paths; configure default figure text

off_range = -100:100; % range of offsets to search for block matches
blksz     = 50 ;    % size of blocks to divide images into
z_range   = -30:30;   % range around best matches to use
auto_crop = true;
% force block size to be even
blksz = blksz + mod(blksz,2);

% load the z stack into a 3d array
stack_names = dir(fullfile(data_dir,subj,'2015-*-*__post_stack_00*_AVERAGE.tif'));
stack_date  = stack_names(1).name(1:10);
stack_path  = fullfile(data_dir, subj, stack_names(1).name);
stack       = load_stack(stack_path,auto_crop);

% divide into blocks (define centers in reference to the stack images)
stack_info  = imfinfo(stack_path);
blk_cntrs   = make_blocks(stack_info(1).Height, stack_info(1).Width, blksz);
nblks       = length(blk_cntrs.x);

% get a list of movies each corresponding to a day of imaging
movie_names = dir(fullfile(data_dir,subj,['2015-*__L01__AVERAGE.tif']));
movie_names = {movie_names(:).name};
ndays       = length(movie_names);

if nargin < 2
    days = 1:ndays;
end
days = days(days<=ndays);
if nargin < 3 
    timepoints = 1;
end

% pre-allocate best z choice array
blks_best_z = nan(nblks,ndays);
for di = days
    
    % pre-allocate correlations array 
    C = nan(nblks, length(off_range), length(off_range), length(stack_info));
    
    % load first frame of this day's movie
    img_path = fullfile(data_dir,subj,movie_names{di});
    img_name = movie_names{di};
    img_date = img_name(1:10);
    
    cmat_dir  = fullfile(data_dir, subj, [img_name(1:end-4) '_corrs']);
    if ~exist(cmat_dir), mkdir(cmat_dir); end 

    % let us know what we're image we're working on
    fprintf('\n\nworking on %s...', img_name);

    for ti = timepoints

        fprintf('\n\tframe %i',ti);
        try
            img = load_and_norm_img(img_path, ti);
        catch
            break
        end
        [imgy imgx] = size(img);

        
        
        % initialize best z with middle of stack
        mid_z = floor(size(stack,3) / 2);
        blks_best_z(:,:) = mid_z;
        

        % get cross correlation for each block
        for bi = 1:nblks

            % create block indices around this block center
            [blk_yix blk_xix] = blk_ctr2ix(blk_cntrs, bi, blksz, imgx, imgy);

            % select block
            block   = img(blk_yix,blk_xix);
            
            % determine the neighborhood of z-values to look at
            z_to_check = blks_best_z(bi,di) + z_range;
            z_to_check = z_to_check(z_to_check > 0 & z_to_check <= size(stack,3));

            % get the cross correlation of this block against the stack image
            [c yoff xoff]   = get_xyz_similarity(stack(:,:,z_to_check), block, blk_cntrs.x(bi), blk_cntrs.y(bi), off_range);
            C(bi,:,:,z_to_check) = c;

            % update best z with the z that has the max correlation
            [~, maxix] = max(c(:));
            [~, ~, zix] = ind2sub(size(c),maxix);
            fprintf('\n\tfavorite z for this block in neighborhood around %i...', blks_best_z(bi, di));
            blks_best_z(bi,di) = z_to_check(zix); 
            if di < ndays
                blks_best_z(bi,di+1) = blks_best_z(bi,di);
            end
            fprintf('s %i', blks_best_z(bi, di));
            rpt = blk_reg_report(img, stack, stack_date, img_date, blk_yix, blk_xix, bi , squeeze(C(bi,:,:,:)), off_range);
            saveas(rpt, fullfile(cmat_dir, sprintf('%s_%s_frame%03i_block%03i_neighbors.png',subj,img_date,ti,bi)));
            
        end
        best_corr_z = squeeze(max(max(C,[],2),[],3));
        all_z = figure(2); clf;
        plot(1:size(best_corr_z,2), best_corr_z);
        
        saveas(all_z, fullfile(cmat_dir, sprintf('%s_%s_frame%03i_summary.png',subj,img_date,ti)));
        % now let's save this correlation matrix
        cmat_name = sprintf('%s_%s_frame%03i.mat',subj, img_date, ti);
        savepath = fullfile(cmat_dir,cmat_name); 
        save(savepath, 'C','best_corr_z','off_range','blk_cntrs', 'blksz')
    end
end


end
