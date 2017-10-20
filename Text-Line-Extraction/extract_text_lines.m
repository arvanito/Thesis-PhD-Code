%%  [img_medial_seams, img_sep_seams] = extract_text_lines(img, smooth, s, sigma, off):
%%
%%  function that runs our proposed method for text line extraction using 
%%  projection profile matching and constrained seam carving.
%%
%%  Input:
%%      img: original image
%%      smooth: smoothing parameter for projection profile matching
%%      sigma: standard deviation of the gaussian filter for energy map computation
%%      s: number of image slices for projection profile matching
%%      off: parameter that corresponds to seam thickness
%%
%%  Output:
%%      img_medial_seams: original image overlaid with medial seams
%%      sep_seams: 2-D matrix where each column contains the coordinates of a separating seam
%%      img_sep_seams: original image overlaid with separating seams
%%
%%  Parameter tuning:
%%    The two main parameters of the function are "smooth" and "s". Both of them are crucial to 
%%    the accuracy of the algorithm. 
%%      - The smoothing is important, because the projection profile 
%%      is very ragged in any real manuscript page, and the amount of smoothing influences the number 
%%      of detected lines. Usually, a very small value is needed in the range [0.01, 0.0001].
%%      - The number of slices is also important and it depends on the image resolution and layout. 
%%      Usually, a slice number of 4 works well in practice. For the dataset of Saabni etal. we use s = 4.
%%      For the manuscript of Aline we use s = 8, since the resolution is much higher.
%%    In general, there is no automatic way to select these parameters.  
%%      - "sigma" is not crucial and small values (<3) usually work well.
%%      - small values of "off" are good for visualization purposes (range [1,5]).
%%
%%
%%  Author: Nikolaos Arvanitopoulos (nick.arvanitopoulos@epfl.ch), 2014
%%

function [img_medial_seams, sep_seams, img_sep_seams, img_final] = extract_text_lines(img, smooth, s, sigma, off)

% convert to grayscale if the original is a color image
if (~ismatrix(img))
    img_gray = rgb2gray(img);
    img_gray = im2double(img_gray); % convert to double
else
    img_gray = im2double(img);  % convert to double
end

% compute edge image using Sobel edge detector
img_edge = edge(img_gray);

% medial seam computation with projection profile matching
fprintf('Computing text line approximations....\n');
[~, medial_seams_inds, L, img_medial_seams] = compute_medial_seams(img, img_edge, smooth, s, off);

% separating seam computation with constrained seam carving
fprintf('Generating minimum energy seams...\n');
[sep_seams, img_sep_seams] = compute_separating_seams(img, img_gray, L, sigma, off);

% original image overlaid with both types of seams
img_final = overlay_medial_seams(img_sep_seams, medial_seams_inds, off);

end