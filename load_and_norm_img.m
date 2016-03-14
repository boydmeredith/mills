function img = load_and_norm_img(img_path,zix)
% img = load_and_norm_img(imgpath)
% read image, subtract its minimum from all pixels, normalize it, and return it
if nargin < 2
    zix = 1;
end
    img = imread(img_path,zix);
    img = double(img - min(img(:)));
    img = normalizeToZeroOne(img);
end
