function stack = load_stack(stack_path)
% stack = load_stack(stack_path)
% 
% load stack into a 3d matrix

tiff_info = imfinfo(stack_path);

nframes = numel(tiff_info);    
stack = zeros(tiff_info(1).Height, tiff_info(1).Width, nframes);      
    for fi = 1:nframes
        this_slice = imread(stack_path,'Index',fi,'Info',tiff_info);
        this_slice = double(this_slice - min(this_slice(:)));
        this_slice = normalizeToZeroOne(this_slice);
        stack(:,:,fi) = this_slice;
        
    end
end

