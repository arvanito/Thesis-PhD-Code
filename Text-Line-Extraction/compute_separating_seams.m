%%  [sep_seams, img_sep_seams] = compute_separating_seams(img, img_gray, L, sigma, off):
%%
%%  function that computes horizontal minimum energy separating seams using constrained 
%%  seam carving. The input to the algorithm is the coordinates of the medial seams 
%%  found in the previous text using the function match_hists.
%%  The seam computation is constrained before two consecutive medial seams.
%%  The procedure is repeated until all pairs of text lines are processed. 
%%
%%  Input:
%%      img: original color image
%%      img_gray: original grayscale image
%%      L: medial seams
%%      sigma: standard deviation of the gaussian filter for energy map computation
%%      off: parameter that corresponds to seam thickness
%%
%%  Output:
%%      sep_seams: coordinates of computed separating seams
%%      img_sep_seams: original color image overlaid with the computed separating seams
%%  
%%  Parameter tuning:
%%    The two main parameters of the function is "sigma" and "off". 
%%      - "sigma" is used for prior blurring 
%%      before computing the gradient energy map. It depends on the image resolution. However, 
%%      it is not a crucial parameter, and small values (<3) work in general well.
%%      - "off" corresponds to seam thickness when the separating seams are presented in a figure.
%%
%%
%%  Author: Nikolaos Arvanitopoulos (nick.arvanitopoulos@epfl.ch), 2014
%%

function [sep_seams, img_sep_seams] = compute_separating_seams(img, img_gray, L, sigma, off)

% blur the image with gaussian filter
img_blur = anigauss_mex(img_gray, sigma);

% compute the energy map
[Ix,Iy] = gradient(img_blur);
energy_map = abs(Ix) + abs(Iy);

% transpose the energy map, because we compute vertical seams
energy_map = energy_map';

n = size(energy_map,1);

% initialize the seam image
img_sep_seams = img;

% number of text lines in the page
l = length(L);

% initialize seam coordinates
sep_seams = zeros(n,l-1);

% copy the energy map for the dynamic programming part 
new_energy = energy_map;

%% apply constrained seam carving for each pair of text lines
for i = 1:l-1
    % upper and lower medial seams
    L_a = L{i};
    L_b = L{i+1};
    
    % length of the medial seam (maybe they do not extend
    % through the whole page width)
    l_a = length(L_a);
    l_b = length(L_b);
    min_l = min(l_a,l_b);   % minimum of the two text line lengths
    
    %% compute minimum energy separating seam using dynamic programming
    for row = 2:min_l
        for col = L_a(row):L_b(row)
            % find previous row's neighbors for the cumulative matrix
            % computation and take care not to overstep the boundaries of 
            % the image
            left = max(col-1,L_a(row-1));
            right = min(col+1,L_b(row-1));

            % find the minimum value of the previous row's neighbors
            minpath = min(new_energy(row-1,left:right));
            
            % take care of discontinuous boundaries in the text line 
            % approximations
            if (isempty(minpath))   % we are in the boundary 
                % many pixels difference, discontinuity
                if (col > left)     
                    new_energy(row,col) = new_energy(row-1,right);
                elseif (col < right)   
                    new_energy(row,col) = new_energy(row-1,left);
                end
            else    % one pixel difference, no discontinuity
                new_energy(row,col) = energy_map(row,col) + minpath;
            end
        end
    end

    %% trace the path backwards
    % find the minimum energy path at the bottom, the index of 
    % the minimum energy is the starting point
    [~, min_index] = min(new_energy(min_l,L_a(min_l):L_b(min_l)));
    min_index = L_a(min_l) + min_index - 1; % correct index in the entire image

    % backtrack through the energy map from bottom to top  
    % to determine the minimum energy seam
    for row = min_l-1:-1:1
        % the column of the minimum energy seam
        j = min_index(1);
        
        % update the seam map
        sep_seams(row+1,i) = j;
        
        % find next row's neighbors
        left = max(j-1,L_a(row));
        right = min(j+1,L_b(row));
        
        % the column of the minimum energy seam
        [~, min_index] = min(new_energy(row,left:right));
        
        % take care of discontinuous boundaries in the text line 
        % approximations
        if (isempty(min_index))     % many pixels difference, discontinuity
            if (j > left)
                min_index = right;
            elseif (j < right)
                min_index = left;
            end
        else    % one pixel difference, no discontinuity
            min_index = min_index + left - 1;   % correct index in the entire image
        end
    end

    % first position of the optimal seam
    sep_seams(1,i) = min_index(1);
    
    %% overlay separating seams on the original image
    img_sep_seams = overlay_separating_seams(img_sep_seams, sep_seams(:,i), off);
    
    % reset the energy values
    new_energy = energy_map; 
end

end