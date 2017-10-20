%%  img_medial_seams = overlay_medial_seams(img, medial_seams_inds, off):
%%
%%  function that overlays the medial_seams on the original 
%%  image, either grayscale or color.
%%
%%  Input:
%%      img: original image, grayscale or color
%%      medial_seams_inds: linear indices of the medial seams
%%      off: parameter that corresponds to seam thickness
%%
%%  Output:
%%      img_medial_seams: original image overlaid with the text line approximations
%%
%%  Parameter tuning:
%%  - "off" corresponds to seam thickness when the separating seams are presented in a figure.
%%
%%
%%  Author: Nikolaos Arvanitopoulos (nick.arvanitopoulos@epfl.ch), 2014
%%

function img_medial_seams = overlay_medial_seams(img, medial_seams_inds, off)

n = size(img,1);
m = size(img,2);

% initialize text lines image
img_medial_seams = img;

% compute (x,y) coordinates of the linear text line approximation
% indices
[I,J] = ind2sub([n,m],medial_seams_inds);

% grayscale image
if (ismatrix(img))
    for ii = 1:length(I)
        % use black color
        if (I(ii) <= off )  % row coordinates very near the first image row
            img_medial_seams(I(ii):I(ii)+off,J(ii)) = 0;
        elseif (I(ii) > n-off)  % row coordinates very near the last image row
            img_medial_seams(I(ii)-off:I(ii),J(ii)) = 0;
        else    % intermediate row coordinates
            img_medial_seams(I(ii)-off:I(ii)+off,J(ii)) = 0;
        end
    end
else    % color image
    for ii = 1:length(I)
        % use blue color
        if (I(ii) <= off )  % row coordinates very near the first image row
            img_medial_seams(I(ii):I(ii)+off,J(ii),1) = 0;
            img_medial_seams(I(ii):I(ii)+off,J(ii),2) = 0;
            img_medial_seams(I(ii):I(ii)+off,J(ii),3) = 255;
        elseif (I(ii) > n-off)  % row coordinates very near the last image row
            img_medial_seams(I(ii)-off:I(ii),J(ii),1) = 0;
            img_medial_seams(I(ii)-off:I(ii),J(ii),2) = 0;
            img_medial_seams(I(ii)-off:I(ii),J(ii),3) = 255;
        else    % intermediate row coordinates
            img_medial_seams(I(ii)-off:I(ii)+off,J(ii),1) = 0;
            img_medial_seams(I(ii)-off:I(ii)+off,J(ii),2) = 0;
            img_medial_seams(I(ii)-off:I(ii)+off,J(ii),3) = 255;
        end
    end
end

end