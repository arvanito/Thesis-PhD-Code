%%  [local_max, medial_seams_inds, L, img_medial_seams] = compute_medial_seams(img, img_edge, smooth, s, off):
%%
%%  function that is based on projection profile matching for medial seam computation.
%%  The image is split into slices and local maxima of the edge image are matched 
%%  between two consecutive slices. After the process of all slices, we obtain 
%%  piece-wise linear approximation of text lines.
%%  The output is a matrix where each element (i,j) contains the column 
%%  coordinate of the i-th seam. j ranges from 1 to the page width.
%%
%%  Input:
%%      img: original image
%%      img_edge: edge image
%%      smooth: smoothing parameter for projection profile matching
%%      s: number of image slices for projection profile matching
%%      off: parameter that corresponds to seam thickness
%%
%%  Output:
%%      local_max: local maxima of the projection profiles of each slice
%%      medial_seams_inds: linear indices of the text line approximations
%%      L: cell array where each cell element contains the row coordinates of the medial seams
%%      img_medial_seams: original image overlaid with the medial seams
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
%%      - small values of "off" are good for visualization purposes (range [1,5]).
%%
%%
%%  Author: Nikolaos Arvanitopoulos (nick.arvanitopoulos@epfl.ch), 2014
%%

function [local_max, medial_seams_inds, L, img_medial_seams] = compute_medial_seams(img, img_edge, smooth, s, off)

% parameter initialization 
[n,m] = size(img_edge);

% binary image that will contain the piece-wise text line approximations
img_bin = zeros(n*m,1);

% width of each image slice
w = floor(m / s);

%% compute horizontal projection profiles of all edge image slices 
%% and find their local maxima
local_max = cell(s,1);   % locations of local maxima
k = 1;  % help counter for the slices
for i = 1:s     % for each slice
    % compute histogram of sum of edges
    horiz_proj = sum(img_edge(:,k:k+w-1),2);
    
    % smooth it with cubic splines
    horiz_proj_smooth = fnval(csaps(1:n,horiz_proj,smooth),1:n);
    
    % find peaks of the profile
    [~, local_max{i}] = findpeaks(horiz_proj_smooth);
    
    % next slice
    k = k + w;
end

% half of slice width <-!!
k = floor(w/2);
medial_seams_inds = [];

%% match local maxima of the projection profiles 
%% between two consecutive image slices
for i = 1:s-1

    % if there is no local max, continue to the next pair
    if (isempty(local_max{i}))
        continue;
    end
    
    % matches from left to right and from right to left
    matches_lr = zeros(length(local_max{i}),1);
    matches_rl = zeros(length(local_max{i+1}),1);
    
    % find matches from left slice to right slice
    for j = 1:length(local_max{i})  % for each left local maximum
        % compute distances from the right slice and sort them
        dists = abs(local_max{i}(j) - local_max{i+1});
        [~,ind_sort] = sort(dists,'ascend');
        
        % the match is the right maximum that is closest to the left one
        matches_lr(j) = ind_sort(1);
    end

    % find matches from right slice to left slice
    for j = 1:length(local_max{i+1})    % for each right local maximum
        % compute distances from the left slice and sort them
        dists = abs(local_max{i+1}(j) - local_max{i});
        [~,ind_sort] = sort(dists,'ascend');
        
        % the match is the left maximum that is closest to the right one
        matches_rl(j) = ind_sort(1);
    end
    
    % match profile maxima that agree from both sides
    % (left to right and right to left)
    for j = 1:length(local_max{i})
        for o = 1:length(local_max{i+1})
            
            % maximum match exists
            if (matches_lr(j) == o && matches_rl(o) == j)
                % calculate rounded row coordinates of lines between 
                % matched maxima
                if (i == 1) % start the first histogram profile from the beginning of the image
                    % rounded row coordinates for the first pair of slices
                    points = round(linspace(local_max{i}(j),local_max{i+1}(o),k));
                    
                    % linear indices on the entire image
                    inds = sub2ind([n,m],points,1:k);
                    medial_seams_inds = [medial_seams_inds, inds];
                elseif (i == s-1) % end the last histogram profile in the end of the image
                    len = length(k:m);
                    
                    % rounded row coordinates for the last pair of slices
                    points = round(linspace(local_max{i}(j),local_max{i+1}(o),len));
                    
                    % linear indices on the entire image
                    inds = sub2ind([n,m],points,k:m);
                    medial_seams_inds = [medial_seams_inds, inds]; %#ok<*AGROW>
                else    % intermediate pairs of slices
                    
                    % rounded row coordinates for an intermediate pair of
                    % slices
                    points = round(linspace(local_max{i}(j),local_max{i+1}(o),w));
                    
                    % linear indices on the entire image
                    inds = sub2ind([n,m],points,k:k+w-1);
                    medial_seams_inds = [medial_seams_inds, inds];
                end
                
                % update the binary image with the text line approximation
                img_bin(inds,:) = 1;
            end
        end
    end
    
    % update the help counter 
    if (i > 1 && i < s-1)
        k = k + w;
    end
end

% reshape the binary image of medial seams
img_bin = reshape(img_bin,[n,m]);

% unique linear indices
medial_seams_inds = unique(medial_seams_inds);

%% post-processing step to remove lines that start from some intermediate 
%% column of the image
% connected component analysis of the binary image containing the 
% piece-wise medial seams
CC = bwconncomp(img_bin);

% number of connected components
num_cc = CC.NumObjects;

img_bin = img_bin(:);

% remove medial seams that do not start from the beginning of
% the image
for c = 1:num_cc
    % column coordinates of the c-th connected component
    inds_cc = CC.PixelIdxList{c};
    [~,J] = ind2sub([n,m],inds_cc);
        
    % if the connected component does not start from the beginning of the 
    % image, remove it from the computed indices
    if (J(1) ~= 1)
        img_bin(inds_cc) = 0;   % update the binary image
        f = ismember(medial_seams_inds,inds_cc);
        medial_seams_inds(f) = [];  % remove the indices
    end
end
img_bin = reshape(img_bin,[n,m]);

%% post-processing step to extend the small lines towards the end 
%% column of the image
% connected component analysis of the binary image containing the 
% piece-wise medial seams
CC = bwconncomp(img_bin);

% number of connected components
num_cc = CC.NumObjects;

% extend medial seams to the end of the image, if possible
for c = 1:num_cc
    % column coordinates of the c-th connected component
    inds_cc = CC.PixelIdxList{c};
    [I,J] = ind2sub([n,m],inds_cc);
    
    % check if this connected component ends in the last column of the 
    % image, if not, extend it towards the end
    if (J(end) ~= m)
        % last column of the CC
        end_col = J(end);
        if (c == 1) % the top connected component
            % create linear indices between the first row and the row 
            % coordinates of the connected component
            inds_ext = floor(linspace(I(end_col),1,m-end_col+1));
            inds_ext = sub2ind([n,m],inds_ext,end_col:m);
            
            % add the extended indices to the full indices array
            medial_seams_inds = [medial_seams_inds, inds_ext];
        elseif (c == num_cc)    % the bottom connected component
            % create linear indices between the last row and the row 
            % coordinates of the connected component
            inds_ext = floor(linspace(I(end_col),n,m-end_col+1));
            inds_ext = sub2ind([n,m],inds_ext,end_col:m);
            
            % add the extended indices to the full indices array
            medial_seams_inds = [medial_seams_inds, inds_ext];
        else    % intermediate connected component
            % extract the lower (previous) CC
            prev_cc = CC.PixelIdxList{c-1};
            [I_p,~] = ind2sub([n,m],prev_cc);
            
            % extract the upper (next) CC
            next_cc = CC.PixelIdxList{c+1};
            [I_n,~] = ind2sub([n,m],next_cc);
            
            % extend the CC in the middle between the row coordinates 
            % of the upper and lower CC's until an intersection is found
            middle = floor((I_p(end)+I_n(end))/2);  % middle end point
            
            % extend towards the end
            inds_ext = floor(linspace(I(end_col),middle,m-end_col+1));
            
            % find intersections of the extended CC to the lower and upper 
            % ones
            [inter_p,ip1,ip2] = intersect(inds_ext,I_p);
            [inter_n,in1,in2] = intersect(inds_ext,I_n);
            if (~isempty(inter_p) && ~isempty(inter_n))   % if both intersections
                if (inds_ext(ip1(1)) < inds_ext(in1(1)))    % first intersection is with the upper CC
                    % extend the approximation until the first intersection
                    t = ip2(1) - end_col;
                    inds_ext = inds_ext(1:t-10);
                    inds_ext = sub2ind([n,m],inds_ext,end_col:end_col+t-1-10);
                else    % first intersection is with the lower CC
                    % extend the approximation until the first intersection
                    t = in2(1) - end_col;
                    inds_ext = inds_ext(1:t-10);
                    inds_ext = sub2ind([n,m],inds_ext,end_col:end_col+t-1-10);
                end
            elseif (~isempty(inter_p))  % intersection with the upper CC
                % extend the approximation until the first intersection
                t = ip2(1) - end_col;
                inds_ext = inds_ext(1:t-10);
                inds_ext = sub2ind([n,m],inds_ext,end_col:end_col+t-1-10);
            elseif (~isempty(inter_n))  % intersection with the lower CC
                % extend the approximation until the first intersection
                t = in2(1) - end_col;
                inds_ext = inds_ext(1:t-10);
                inds_ext = sub2ind([n,m],inds_ext,end_col:end_col+t-1-10);
            else    % if no intersection, just convert them to linear indices
                inds_ext = sub2ind([n,m],inds_ext,end_col:m);
            end

            % add the extended indices to the full indices array
            medial_seams_inds = [medial_seams_inds, inds_ext];
        end
    end
end

% binary image with only the medial seams
img_bin = zeros(n*m,1);
img_bin(medial_seams_inds) = 1;
img_bin = reshape(img_bin,n,m);

% connected component analysis of the final binary image containing 
% the medial seams
CC = bwconncomp(img_bin);
num_cc = CC.NumObjects;

% smooth the piece-wise medial seams with a spline 
% and save them 
L = cell(num_cc,1);

% for each connected component
for c = 1:num_cc    
    inds_cc = CC.PixelIdxList{c};
    
    % discard small components, here threshold!!!!!
    if (length(inds_cc) < 10)   
        continue;
    end
    
    % calculate (x,y) coordinates of connected component
    [I,J] = ind2sub([n,m],inds_cc);
    
    % smooth the component with a spline
    sp = spaps(J,I,0.1);
    
    % evaluate the spline and round the output
    L{c} = round(fnval(sp,J));
end

% overlay the medial seams on the original image
img_medial_seams = overlay_medial_seams(img, medial_seams_inds, off);

end