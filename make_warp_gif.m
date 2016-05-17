function [im] = make_warp_gif(cmat_dir, cmat_names, gifname)
% make_warp_gif(cmat_dir, cmat_names)
% 
% return movie visualizing a manifold in 3d space describing the 
% deformation of a brain during imaging by plotting
% the z value that maximizes cross correlation for 
% an image block at the relevant xy position
%
% cmat_dir      a root directory in which to look 
%               for cmat files (correlation matrices 
%               computed previously) and save the gif
% cmat_names    names of those cmat files with the  
%               cmat_dir

ndays   = length(cmat_names);
gifname = fullfile(cmat_dir, [gifname '.gif']);
b = listfonts;
set(0,'defaulttextfontname',b{1});
for di = 1:ndays
    cmat_path = fullfile(cmat_dir,cmat_names{di});
    if exist(cmat_path,'file')
        load(cmat_path, 'best_corr_z')
    end

    [~,maxZ] = max(best_corr_z,[],2);
    grid_size = sqrt(size(maxZ,1));
    assert(0==mod(grid_size,1));
    maxZ = reshape(maxZ, grid_size, grid_size);
    maxZ = maxZ([2:(grid_size-1)],[2:(grid_size-1)]);    
    h = figure(1);clf;surf(maxZ)
    set(h,'color','w');
    zlim([1 51]);
    colormap redbluecmap;
    grid off;
    title(strrep(cmat_names{di},'_',' '),'FontName','Arial'); 
    colorbar;
    caxis([24 42 ]);
    set(gca,'zdir','reverse');
    f.cdata = hardcopy(h, '-Dzbuffer', '-r0');
    if di == 1
        [im, map] = rgb2ind(f.cdata,256,'nodither');
        im(1,1,1,ndays) = 0;
    else
        im(:,:,1,di) = rgb2ind(f.cdata,map,'nodither');
    end
end

if ~isempty(gifname)
    imwrite(im,map,gifname,'DelayTime',.2,'LoopCount',Inf);
end