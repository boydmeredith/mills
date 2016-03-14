function img = load_and_norm_img(img_path)
% img = load_and_norm_img(imgpath)
% read image, subtract its minimum from all pixels, normalize it, and return it

    img = imread(img_path);
    img = double(img - min(img(:)));
    img = normalizeToZeroOne(img);
end
