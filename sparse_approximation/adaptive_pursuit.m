function [in] = adaptive_pursuit(in)

% Adaptive Pursuit - die eierlegende Wollmilchsau
%
% Input:
%   dico: dictionary matrix 
%   Y: matrix: columns are signals we want to approximate
%   maxit: maximal number of iterations 
%                   (default max. sensible sparsity level floor(d/log(K)))
%   precision: smallest coefficient size (default - 10^machine precision)
%
% Output:
%   X: matrix: columns are coefficients of the s-sparse approximations
%               of the signals in the dictionary
%   normres: relative residual norms
%   I: Indices of the supports of X
%   iter: number of iterations

if (~strcmp(in.SPARSITY_MODE, 'ADAPTIVE'))
    error('This implementation of adaptive pursuit sparse recovery needs adaptive sparsity!');
end

if (~isfield(in, 'SPARSITY_NOISE_RISK'))
    in.SPARSITY_NOISE_RISK = 1/4;
    warning(['SPARSITY_NOISE_RISK is set to ' num2str(in.SPARSITY_NOISE_RISK)]);
end

in.X = zeros(in.K, in.N);
iter = zeros(1, in.N);
sparsities = zeros(1, in.N);
normres = zeros(1, in.N);
I = zeros(in.K, in.N);

maxit = floor(in.D/log(in.K)); %%% maximum number of iterations

precision = 10^(-8); %%% could be a parameter

th = sqrt(2*(log(2*in.K)-log(in.SPARSITY_NOISE_RISK))/in.D);

number_iterations_done = zeros(in.N, 1);

for n = 1:in.N 
    
    y = in.Y(:, n);
    
    coeff = zeros(in.K, 1);
    app = zeros(in.D, 1);
    res = y;
    
    %supp_ancient = [];
    supp_old = -1;
    supp_new = [];
   
    it = 0;
    
    while length(setxor(supp_old,supp_new)) > 0 && it < maxit
        supp_old = supp_new;
        
        ip = coeff + in.DICT'*res; 
%         supp_new = find(abs(ip) > th*norm(res) + precision);

%         app = in.DICT(:, supp_new)*coeff;
%         res = in.Y(:,n) - app;
        enres = res'*res;
        enapp = app'*app;
    
        supp_new = find(abs(ip) > sqrt(enapp/in.D+enres*2*log(in.K/in.SPARSITY_NOISE_RISK)/in.D));
        
        suppsize = length(supp_new);
        if (suppsize > in.S_MAX) || (suppsize < in.S_MIN)
            [~, sorted_indices] = sort(abs(ip), 'descend');
            if (suppsize > in.S_MAX)
                supp_new = sorted_indices(1:in.S_MAX);
            end

            if (suppsize < in.S_MIN)
                supp_new = sorted_indices(1:in.S_MIN);
            end
        end
        
        coeff = zeros(in.K, 1);
        coeff(supp_new) = in.DICT(:, supp_new)\y;
        app = in.DICT*coeff;
        res = y - app;
        
        it = it+1;
    end
    
    number_iterations_done(n) = it;
    
    in.X(:, n) = coeff;
    sparsities(n) = length(supp_new);
    normres(n) = norm(res)/norm(y);
    I(supp_old, n) = 1;
    iter(n) = it;
end

% mean(number_iterations_done)
mean(sparsities)
