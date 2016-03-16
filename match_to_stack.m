function best_z = match_to_stack(stack, slice_img, varargin)
setup; 
pairs = {'resize_by'   1;...
         'z_range'     [0 0];...
         'eq_hist'     1;...
         'max_d'       25;...
}; parseargs(varargin, pairs);

n_z       = size(stack, 3);
corr      = nan(1,n_z); 
d         = nan(1,n_z);
z_check   = (1 + z_range(1)):(n_z - z_range(2));

slice_hist = imhist(slice_img);
slice_img = imresize(slice_img,resize_by);


for zi = z_check
    zi_best_corr = -Inf;
    zi_chosen_d  = nan;
    
    temp = stack(:,:,zi);
    if eq_hist, temp = histeq(temp, slice_hist); end
    temp = imresize(temp, resize_by);

    [this_max_c this_d] = window_normxcorr2(temp,slice_img,max_d);
        
    if this_max_c > zi_best_corr && this_d.d < max_d
        corr(zi)    = this_max_c;
        d(zi)       = this_d.d;
    end
corr(isnan(corr)) = -Inf;
[~, sorted_z] = sort(corr,'descend');
corr(isinf(corr)) = nan;
best_z = sorted_z(1);
end



end