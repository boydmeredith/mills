cmat_dir = '/Volumes/tank/jtb3/Data/J114/2015-09-25__L01__AVERAGE_corrs/';
cmat_names = dir(fullfile(cmat_dir,'*frame*mat'));
cmat_names = {cmat_names.name};
blocknos = [17 25 33 73 86];
nblocks = length(blocknos);
nframes   = length(cmat_names);
for bi = 1:length(blocknos)
    blockno = blocknos(bi);
    z_v_t =  zeros(26, nframes);
    for di = 1:nframes
        cmat_path = fullfile(cmat_dir,cmat_names{di});
        fprintf('loading frame %i',di);
        if exist(cmat_path,'file')
            load(cmat_path, 'best_corr_z')
        end
        
        z_v_t(:,di) = squeeze(best_corr_z(blockno,:));
    end
    
    [~, maxix] = max(z_v_t);
    freq = zeros(1,26); 
    for i = 1:26
        freq(i) = sum(maxix==i)/size(maxix,2);
    end

    subplot(nblocks,4,[((bi-1)*4+1):((bi-1)*4+3)]); 
    imagesc(z_v_t); hold on;
    plot(maxix,'k'); 
    set(gca,'YDir','normal'); 
    ylabel('z slice #'); xlabel('movie frame #'); 
    title(sprintf('J114 09-25 block %i',blockno)); 
    c  = colorbar('WestOutside'); ylabel(c,'max correlation');

    set(gca,'zdir','reverse');

    subplot(nblocks,4,4*bi); 
    bar(gca, freq, 'FaceColor',[0,0,0]);
    xlim([1 26]); 
    set(gca, 'view', [90 -90]); box off
    %set(gca,'xdir','reverse');
end
