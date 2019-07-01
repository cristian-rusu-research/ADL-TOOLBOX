function [in] = thresholding_adaptive(in)

if (~strcmp(in.SPARSITY_MODE, 'ADAPTIVE'))
    error('This implementation of thresholding sparse recovery needs adaptive sparsity!');
end

if (~isfield(in, 'SPARSITY_NOISE_RISK'))
    in.SPARSITY_NOISE_RISK = 1/4;
    warning(['SPARSITY_NOISE_RISK is set to ' num2str(in.SPARSITY_NOISE_RISK)]);
end

if (~isfield(in, 'S_ESTIMATED'))
    in.S_ESTIMATED = max(in.S_MIN, 1);
    warning(['S_ESTIMATED is set to S_MIN value of ' num2str(in.S_ESTIMATED)]);
end

ip = in.DICT'*in.Y;
absip = abs(ip);
[~, I] = sort(absip, 1, 'descend');
G = in.DICT'*in.DICT;

in.SPARSITY = in.S_ESTIMATED*ones(1, in.N);
Sguess = 0;
for n=1:in.N
    In = I(1:in.S_ESTIMATED,n);
    coeff = G(In,In)\ip(In,n);
    
    %%% if the line above creates a warning (instable)
    %%% recalculate stably but more expensively via pinv
    [~, msgid] = lastwarn;
    if ~isempty(msgid)
        if strcmp(msgid, 'MATLAB:nearlySingularMatrix') || ...
            strcmp(msgid, 'MATLAB:singularMatrix')
            coeff = pinv(G(In,In))*ip(In,n);
        end
    end
    
%%%%%% sequential estimation of sparsity level 
% %     app = in.DICT(:,In)*coeff;
% %     res = in.Y(:,n) - app;
% %     enres = res'*res;
% %     enapp = app'*app;
% %     
% %     resip = in.DICT(:,I(in.S_ESTIMATED+1:in.K, n))'*res;
% %     Sguess = Sguess + length(find(abs([coeff;resip]) > sqrt(enapp/in.D+enres*2*log(in.K/in.SPARSITY_NOISE_RISK)/in.D)));
      
    in.X(:, n) = 0;
    in.X(In, n) = coeff;
end

%%% sequential estimation of sparsity level
%%% Sguess = Sguess/in.N

%%% fast estimation of sparsity level
R = in.Y - in.DICT*in.X;
% the residual norms
residual_energy = sum(R.^2); 
app_energy = sum((in.Y).^2)-residual_energy;
sig_threshold = sqrt((app_energy+residual_energy*2*log(in.K/in.SPARSITY_NOISE_RISK))/in.D);
% inner products with residual + coefficients
Z = in.X + in.DICT'*R;   

Sguess = sum(sum(ceil(abs(Z) - ones(in.K,1)*sig_threshold)))/in.N;

% compute the sparsity level for the next iteration
in.S_ESTIMATED = in.S_ESTIMATED + sign((round(Sguess)) - in.S_ESTIMATED);
in.S_ESTIMATED = min([in.S_MAX, in.S_ESTIMATED]);
in.S_ESTIMATED = max([in.S_ESTIMATED, in.S_MIN]);

end
