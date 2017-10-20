%%  run_text_line_extraction(dir_in, dir_out, file_ext):
%%
%%  function that runs our text line extraction algorithm on 
%%  all images of an input directory with a specific file extension.
%%  All standard image extensions are valid, such as tiff, jpg, png.
%%
%%  Example of use:
%%
%%  dir_in = '/Users/Nikolaos/Documents/MATLAB/input_images/';
%%  dir_out = '/Users/Nikolaos/Documents/MATLAB/results/';
%%  file_ext = 'jpg';
%%
%%  run_text_line_extraction(dir_in, dir_out, file_ext);
%%
%%  Author: Nikolaos Arvanitopoulos (nick.arvanitopoulos@epfl.ch), 2014
%%

function run_text_line_extraction(dir_in, dir_out, file_ext)

% list of files in the directory name with the input file extension
listing = dir(strcat(dir_in,'*.',file_ext));
file_names = {listing.name};

% number of pages in the directory with this file extension
num_pages = length(file_names);

% parameters initialization
smooth = 0.0003;    % smoothing for projection profile matching
s = 8;              % number of image slices for projection profile matching
sigma = 3;          % standard deviation for gaussian smoothing
off = 5;            % offset for seam overlay, corresponds to seam thickness

% process all pages in the directory
for i = 1:num_pages
    fprintf('Processing page No: %d\n', i);
    
    % load the image from the directory
    img = imread(strcat(dir_in,file_names{i}));
    
    % call the text line extraction code
    [img_text_lines, sep_seams, img_seams, img_final] = extract_text_lines(img, smooth, s, sigma, off); %#ok<ASGLU>
    
    % remove the file extension from the file name
    name = strrep(file_names{i},strcat('.',file_ext),'');
    
    % create a directory with the file name without the extension
    complete_dir_out = strcat(dir_out,name);
    
    %% create standard extension names for each file image and save them
    
    % save medial seam images
    text_lines_name = strcat(complete_dir_out,'_text_lines','.',file_ext);
    imwrite(img_text_lines,text_lines_name);
    
    % save separating seams in a .mat file
    sep_seams_name = strcat(complete_dir_out,'_seams.mat');
    save(sep_seams_name,'sep_seams');
    
    % save separating seam images
    seams_name = strcat(complete_dir_out,'_seams','.',file_ext);
    imwrite(img_seams,seams_name);
    
    % save images with both seams overlaid
    final_name = strcat(complete_dir_out,'_final','.',file_ext);
    imwrite(img_final,final_name);
end