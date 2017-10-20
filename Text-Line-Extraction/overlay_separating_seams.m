%%  img_sep_seams = overlay_separating_seams(img, sep_seams, min_l, off):
%%
%%  function that overlays computed seams on the original image. It 
%%  takes into account if it is color or grayscale. An offset parameter 
%%  controls the width of the presented seam.
%%
%%  Input:
%%      img: original image (grayscale or color)
%%      sep_seams: seam coordinates, each column contains one seam
%%      off: parameter that corresponds to seam thickness
%%
%%  Output:
%%      img_sep_seams: original image overlaid with the computed seams
%%
%%  Parameter tuning:
%%      - "off" corresponds to seam thickness when the separating seams are presented in a figure.
%%
%%
%%  Author: Nikolaos Arvanitopoulos (nick.arvanitopoulos@epfl.ch), 2014
%%

function img_sep_seams = overlay_separating_seams(img, sep_seams, off)

% initialize seam image
img_sep_seams = img;

% number of seams
l = size(sep_seams,2);
le = size(sep_seams,1);

% grayscale image
if (ismatrix(img))
    % for each seam
    for j = 1:l  
        if (sum(sep_seams(:,j)) > 0)
            for i = 1+off:le-off
                % use black color for the seam
                if (sep_seams(i,j) > 0)
                    img_sep_seams(sep_seams(i,j),i-off:i+off) = 0; 
                end
            end
        end
    end
else    % color image
    % for each seam
    for j = 1:l     
        if (sum(sep_seams(:,j)) > 0)
            for i = 1+off:le-off
                % use red color 
                if (sep_seams(i,j) > 0)
                    img_sep_seams(sep_seams(i,j),i-off:i+off,1) = 255; 
                    img_sep_seams(sep_seams(i,j),i-off:i+off,2) = 0;
                    img_sep_seams(sep_seams(i,j),i-off:i+off,3) = 0;
                end
            end
        end
    end
end


end