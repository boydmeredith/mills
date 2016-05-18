function [] = display_best_z(subj,datestring,z_range,block_range)
%
% show best match of

%%

setup;
if nargin<2
    z_range = -1:1;
end
set(0,'defaultfigurevisible','off');
set(0,'defaultaxesfontsize',12);
stack_names = dir(fullfile(data_dir,subj,'2015-*-*__post_stack_00*_AVERAGE.tif'));
stack_date  = stack_names(1).name(1:10);
stack_path  = fullfile(data_dir, subj, stack_names(1).name);
stack       = load_stack(stack_path);
[stacky stackx stackz] = size(stack);

stack_info  = imfinfo(stack_path);

off_in_cplot = -30:30;

movie_names = dir(fullfile(data_dir,subj,['2015-*' datestring '*__L01__AVERAGE.tif']));
movie_names = {movie_names(:).name};
ndays       = length(movie_names);


for di = 1:ndays
    % load first frame of this day's movie
    img_path = fullfile(data_dir,subj,movie_names{di});
    img_name = movie_names{di};
    img_date = img_name(1:10);
    img = load_and_norm_img(img_path);
    ih = imhist(img);
    [imgy imgx] = size(img);
    
    corr_report_dir = [img_path(1:end-4) '_corr_reports'];
    if ~exist(corr_report_dir,'dir')
        mkdir(corr_report_dir);
    end
    
    
    % load correlations for this day's movie
    corrs_path = [img_path(1:end-4) '_correlations.mat'];
    
    if ~exist(corrs_path)
        continue;
    end
    load(corrs_path);
    nblks = length(blk_cntrs.x);
    
    cplot_ix = find(ismember(off_range,off_in_cplot));
    
    
    for blockno = block_range
        this_block_C    = squeeze(C(blockno,:,:,:));
        % get the block from the stack that corresponds to max correlation for
        % this block
        [minC]          = min(this_block_C(:));
        [maxC, maxix]   = max(this_block_C(:));
        
        [i, j, best_z]  = ind2sub(size(this_block_C), maxix);
        
        xoff = off_range(j);
        yoff = off_range(i);
        
        img_ctr_x = blk_cntrs.x(blockno) ;
        img_ctr_y = blk_cntrs.y(blockno) ;
        
        [img_block_yix img_block_xix] = blk_ctr2ix(blk_cntrs, blockno, blksz,...
            imgx, imgy);
        blk_cntrs_off = blk_cntrs;
        blk_cntrs_off.x(blockno) = blk_cntrs.x(blockno) + xoff;
        blk_cntrs_off.y(blockno) = blk_cntrs.y(blockno) + yoff;
        [stack_block_yix stack_block_xix] = blk_ctr2ix(blk_cntrs_off, blockno,...
            blksz, stackx, stacky);
        
        
        
        % now need to grab
        figure(1); clf
        set(figure(1), 'visible','off','Position', [130    15   1500   1000], 'PaperPositionMode','auto')
        [~,hMat] = makeSubplots(1,4,length(z_range),.01,.1,[.02 0 .97 .98]);
        %     axes(hMat(1,1)); title('stack image','FontSize',15); hold on;
        %     axes(hMat(1,2)); title('movie'); hold on;
        %     axes(hMat(1,3)); title('2d xcorr'); hold on;
        annotation('textbox', [0 0.9 1 0.1], ...
            'String', sprintf('stack date: %s;    movie date: %s;     block: %i',stack_date,img_date, blockno), ...
            'EdgeColor', 'none', ...
            'HorizontalAlignment', 'center', 'FontSize',15)
        
        for zi = 1:length(z_range)
            z = z_range(zi) + best_z;
            
            
            if z > 0 & z <= size(stack,3)
                axes(hMat(zi,1));
                colormap(hMat(zi,2), bone);
                stack_block = squeeze(stack(stack_block_yix,stack_block_xix,z));
                stack_block = histeq(stack_block,ih);
                imagesc(stack_block);
                hold on;
                
                set(gca, 'visible', 'off');
                if z == best_z
                    ylabel(sprintf('z = %i (best match)',z));
                else
                    ylabel(sprintf('z = %i',z));
                end
                
                
                axes(hMat(zi,2));
                colormap(hMat(zi,1), bone);
                img_block = img(img_block_yix, img_block_xix);
                imagesc(img_block);
                hold on;
                
                set(gca, 'visible', 'off');
                
                axes(hMat(zi,3));
                set(gca, 'visible', 'off');
                if isequal(size(stack_block),size(img_block))
                    imagesc(cat(3,stack_block,img_block,img_block));
                end

                if ~isnan(this_block_C(1,1,z));
                    axes(hMat(zi,4));
                    colormap(hMat(zi,4), hot(20));
                    contourf(this_block_C(cplot_ix,cplot_ix,z));
                    caxis([minC maxC]);
                    hold on;
                    %scatter(j, i, 'ok', 'filled');
                    colorbar
                end
                if zi == length(z_range)
                    set(gca,'XTick', 1:20:length(cplot_ix),...
                        'XTickLabel', off_in_cplot(1:20:length(cplot_ix)));
                else
                    set(gca,'XTick',[])
                end
                set(gca,'YTick', 1:20:length(cplot_ix),...
                        'YTickLabel', off_in_cplot(1:20:length(cplot_ix)));
                
                
                
                
            end
            
        end
                saveas(figure(1),fullfile(corr_report_dir,sprintf('block%02i.png',blockno)))

        best_corr_z = squeeze(max(max(C,[],2),[],3));
        
    end
    
    figure(2); clf; 

    plot(best_corr_z, 1:length(best_corr_z));


    
    keyboard;
end
