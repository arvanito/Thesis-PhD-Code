%   
%   Code for our Laplacian l0 approach for single image reflection removal. 
%   It is built on top of the method described in the following paper 
%   [1] "Image Smoothing via L0 Gradient Minimization", Li Xu, Cewu Lu, Yi Xu, Jiaya Jia, ACM Transactions on Graphics, 
%   (SIGGRAPH Asia 2011), 2011. 

function S = Laplacian_l0(Im, lambda, kappa)
% Laplacian_l0 - Single-Image Reflection removal based on Laplacian l0 Minimization
%   S = Laplacian_l0(Im, lambda, kappa) performs Laplacian-based l0 optimization of input
%   image Im, with regularization parameter lambda and rate kappa.
%
%   Parameters: 
%   Im     : Input uint8 image, both grayscale and color images are acceptable.
%   lambda : Smoothing parameter controlling the degree of smooth. (See [1]) 
%            Typically it is within the range [1e-3, 1e-1], 2e-2 by default.
%   kappa  : Parameter that controls the rate. (See [1])
%            Small kappa results in more iteratioins and with sharper edges.   
%            We select kappa in (1, 2].    
%            kappa = 2 is suggested for natural images.  
%
%   Example of Use
%   ==========
%   Im  = imread('pflower.jpg');
%   S  = Laplacian_l0(Im); % Default Parameters (lambda = 2e-2, kappa = 2)
%   figure, imshow(Im), figure, imshow(S);


% if kappa does not exist, default is 2
if (~exist('kappa','var'))
    kappa = 2.0;
end

% if lambda does not exist, 0.002 is default
if ~exist('lambda','var')
    lambda = 2e-3;
end

% make the image double to interval [0,1]
S = im2double(Im);
I = S;

% RGB or grayscale
D = size(Im,3);

% maximum value of penalizer for the augmented L2 term
betamax = 1e5;

% auxiliary variables
beta = 2*lambda;
it = 0;

% default parameters for adam
alpha = 0.001;
beta_1 = 0.9;
beta_2 = 0.999;
epsilon = 10^(-8);
adam_it = 100;

% run the algorithm
while beta < betamax
    
    % increment the number of iterations
    it = it + 1;
    
    %% H-V subproblem
    % gradients w.r.t x and y direction
    h = [diff(S,1,2), S(:,1,:) - S(:,end,:)];
    v = [diff(S,1,1); S(1,:,:) - S(end,:,:)];
    
    % solve for (H,V) (Eq.(12) of the paper)
    if (D == 1)   % grayscale
        t = ((h.^2 + v.^2) < lambda/beta);
    else        % rgb
        t = (sum((h.^2 + v.^2),3) < lambda/beta);
        t = repmat(t,[1,1,D]);
    end
    h(t)=0; 
    v(t)=0;
    
    
    %% S subproblem
    % handle for the gradient of the objective function 
    % wrt S
    fg = @(S) gradientSGL(S, I, h, v, beta);
    
    % initialize auxiliary variables for Adam
    mt = 0;
    vt = 0;
    
    % run Adam for a fixed number of iterations
    for t = 1:adam_it
        grad = fg(S);
        mt = beta_1 * mt + (1 - beta_1) * grad;
        vt = beta_2 * vt + (1 - beta_2) * (grad .* grad);
        mth = mt / (1 - beta_1^t);
        vth = vt / (1 - beta_2^t);
        S = S - alpha * mth ./ (sqrt(vth) + epsilon);
    end
    
    % increate the penalty of the augmented L2 term
    beta = beta*kappa;
    fprintf('.');
end

fprintf('\n');
end


%% gradient of the objective function w.r.t. S
function grad = gradientSGL(S, I, h, v, beta)

% finite difference operators in x and y
nabla_S_x = [diff(S,1,2), S(:,1,:) - S(:,end,:)];
nabla_S_y = [diff(S,1,1); S(1,:,:) - S(end,:,:)];

% differences between nablas and h, v
nabla_S_x_h = nabla_S_x - h;
nabla_S_y_v = nabla_S_y - v;

% transpose finite difference operators in x and y
nabla_t_S_x_h = [nabla_S_x_h(:,end,:) - nabla_S_x_h(:,1,:), -diff(nabla_S_x_h,1,2)];
nabla_t_S_y_v = [nabla_S_y_v(end,:,:) - nabla_S_y_v(1,:,:); -diff(nabla_S_y_v,1,1)];

% laplacian and transpose
f3 = [0,1,0;1,-4,1;0,1,0];
S_l = zeros(size(S));
I_l = zeros(size(S));
S_l_t = zeros(size(S));
I_l_t = zeros(size(S));
if (ndims(S) == 3)
    for d = 1:3
        S_l(:,:,d) = conv2(S(:,:,d),f3,'same');
        I_l(:,:,d) = conv2(I(:,:,d),f3,'same');
        S_l_t(:,:,d) = conv2(S_l(:,:,d),f3,'same');
        I_l_t(:,:,d) = conv2(I_l(:,:,d),f3,'same');
    end
else
    S_l = conv2(S,f3,'same');
    I_l = conv2(I,f3,'same');
    S_l_t = conv2(S_l,f3,'same');
    I_l_t = conv2(I_l,f3,'same');
end

% final gradient
grad = 2 * (S_l_t - I_l_t) + 2 * beta * (nabla_t_S_x_h + nabla_t_S_y_v);

end

