function stack = load_stack(stack_path, auto_crop)
% stack = load_stack(stack_path)
% 
% load stack into a 3d matrix

if nargin < 2
    auto_crop = false;
end
min_nan_row = Inf;
min_nan_col = Inf;
tiff_info = imfinfo(stack_path);

nframes = numel(tiff_info);    
stack = zeros(tiff_info(1).Height, tiff_info(1).Width, nframes);      
for fi = 1:nframes
    this_slice = imread(stack_path,'Index',fi,'Info',tiff_info);
    this_slice = double(this_slice - min(this_slice(:)));
    this_slice = normalizeToZeroOne(this_slice);
    if auto_crop,
        [y x] = size(this_slice);
        empty_rows = ismember(this_slice, zeros(1,x),'rows');
        empty_cols = ismember(this_slice', zeros(1,y),'rows');
        if min(find(empty_rows)) < min_nan_row
            min_nan_row = min(find(empty_rows));
        end
        if min(find(empty_cols)) < min_nan_col
            min_nan_col = min(find(empty_cols));
        end
        %this_slice(empty_rows,:) = nan;
        %this_slice(:,empty_cols) = nan;

    end

    stack(:,:,fi) = this_slice;
end
if auto_crop
    stack = stack(1:min_nan_row,1:min_nan_col,:);
end
end
