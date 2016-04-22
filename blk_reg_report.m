function [h] = blk_reg_report(img, stack, stack_date,  img_date, img_blk_yix, img_blk_xix, blkno, blk_C, off_range)
%
off_in_cplot = -50:50;
z_range = -3:3;

cplot_ix = find(ismember(off_range,off_in_cplot));
ih = imhist(img);
% get the blk from the stack that corresponds to max correlation for
% this blk
[minC]          = nanmin(blk_C(:));
[maxC, maxix]   = nanmax(blk_C(:));
[i, j, best_z]  = ind2sub(size(blk_C), maxix);

xoff = off_range(j);
yoff = off_range(i);

[stacky stackx stackz] = size(stack);

stack_blk_yix  = img_blk_yix + yoff; 
stack_blk_xix  = img_blk_xix + xoff;
stack_blk_yix  = stack_blk_yix(stack_blk_yix <= stacky & stack_blk_yix > 0);
stack_blk_xix  = stack_blk_xix(stack_blk_xix <= stackx & stack_blk_xix > 0);


h = figure(1); clf
set(figure(1), 'visible','off','Position', [130    15   1500   1000], 'PaperPositionMode','auto')
[~,hMat] = makeSubplots(1,4,length(z_range),.01,.1,[.02 0 .97 .98]);
annotation('textbox', [0 0.9 1 0.1], ...
    'String', sprintf('stack date: %s;    movie date: %s;     blk: %i',stack_date,img_date, blkno), ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center', 'FontSize',15)


for zi = 1:length(z_range)
    z = z_range(zi) + best_z;
    
    
    if z <= 0 || z > size(stack,3)
        continue
    end

    axes(hMat(zi,1));
    colormap(hMat(zi,2), bone);
    stack_blk = squeeze(stack(stack_blk_yix,stack_blk_xix,z));
    stack_blk = histeq(stack_blk,ih);
    imagesc(stack_blk); axis image;
    hold on;
    freezeColors;
    
    set(gca, 'YTick', [], 'XTick',[]);
    hold on;
    if z == best_z
        ylabel(sprintf('z = %i (best match)',z));
    else
        ylabel(sprintf('z = %i',z));
    end
    
    axes(hMat(zi,2));
    colormap(hMat(zi,1), bone);
    img_blk = img(img_blk_yix, img_blk_xix);
    imagesc(img_blk);axis image;
    hold on;
    freezeColors;

    set(gca, 'YTick', [], 'XTick',[]);
    
    axes(hMat(zi,3));
    if isequal(size(stack_blk),size(img_blk))
        imagesc(cat(3,stack_blk,img_blk,img_blk));

        colormap(hMat(zi,3),bone);
    end
    axis image;
    set(gca, 'YTick', [], 'XTick',[]);

    if ~isnan(blk_C(1,1,z));
        axes(hMat(zi,4));
        hold on;
        imagesc(blk_C(cplot_ix,cplot_ix,z));
        colormap redbluecmap
        colorbar
        axis image; 
        freezeColors;
        caxis([minC maxC]);
        %scatter(j, i, 'ok', 'filled');
    end

    if zi == length(z_range)
        set(gca,'XTick', 1:20:length(cplot_ix),...
            'XTickLabel', off_in_cplot(1:20:length(cplot_ix)),'fontsize',10);
    else
        set(gca,'XTick',[])
    end
    set(gca,'YTick', 1:20:length(cplot_ix),...
            'YTickLabel', off_in_cplot(1:20:length(cplot_ix)),'fontsize',10);
end


end
