clc;
clear all;
close all;

% add the path
addpath('method');

% define the input and output paths
inputFolder = '';
outputRoot = '';

% get the images list
imageFiles = dir(fullfile(inputFolder, '*.png'));

%timer
totalTic = tic;

for k = 1:length(imageFiles)
    imgname = imageFiles(k).name;
    imgPath = fullfile(inputFolder, imgname);
    [~, name, ~] = fileparts(imgname);

    img = im2double(imread(imgPath));

    omega = 0.8;
    [J, ~, ~] = physical_enhancement(img, omega);

    output = visibility_enhancement(im2uint8(J));

    result = color_enhancement(J);

    img1 = double(output);
    img2 = double(result);
    imgs_rgb = zeros([size(img1), 2], 'double');
    imgs_rgb(:,:,:,1) = img1;
    imgs_rgb(:,:,:,2) = img2;
    imgs_rgb = imgs_rgb / 255.0;
    
    imgs_lum = rgb2lum(imgs_rgb);
    w1 = get_weight1(imgs_lum);
    w2 = get_weight2(imgs_lum);
    
    w_raw = (w1.^1) .* (w2.^1);
    w = refine_weight(w_raw);
    
    lev = 7;
    img_result = fusion_pyramid(imgs_rgb, w, lev);
    imwrite(img_result, fullfile(outputRoot, [name '.png']));

    clear J output result img1 img2 imgs_rgb imgs_lum w1 w2 w img_result average_result
end

totalTime = toc(totalTic);

