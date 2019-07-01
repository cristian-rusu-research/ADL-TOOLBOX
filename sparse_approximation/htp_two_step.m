function [in] = htp_two_step(in)

if (~strcmp(in.SPARSITY_MODE, 'ADAPTIVE'))
    error('This implementation of 2-step HTP sparse recovery needs adaptive sparsity!');
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
G = in.DICT'*in.DICT;

threshold = sqrt(2*(log(2*in.K)-log(in.SPARSITY_NOISE_RISK))/in.D);   

% update the X, not sparse anymore
Z = in.X + ip - G*in.X;
[~, I] = sort(abs(Z), 1, 'descend');

%% first step, compute best coeff on support of largest S entries
in.SPARSITY = in.S_ESTIMATED*ones(1, in.N);
warning off;
for n = 1:in.N
    In = I(1:in.S_ESTIMATED,n);
    coeff = G(In,In)\ip(In,n);
    
    %%% if the line above creates a warning (instable)
    %%% recalculate stably but more expensively via pinv
    [~, msgid] = lastwarn;
    if ~isempty(msgid)
        if strcmp(msgid, 'MATLAB:nearlySingularMatrix') || ...
            strcmp(msgid, 'MATLAB:singularMatrix')
            coeff = pinv(in.G(In,In))*in.IP(In,n);
        end
    end
    
    in.X(:, n) = 0;
    in.X(In, n) = coeff;
end

if in.K > 1
    %% second step, check which 
    % the residual
    R = in.Y - in.DICT*in.X;
    % the residual norms
    residual_norms = sqrt(sum(R.^2)); 
    % inner products with residual + coefficients
    Z = in.X + in.DICT'*R;   
    %Z = in.X + ip - G*in.X;
    absP = abs(Z);

    [~, I] = sort(absP, 1, 'descend');

    for n = 1:in.N
        In = find(absP(:, n) >= threshold*residual_norms(n));
        in.SPARSITY(n) = length(In);
        if (length(In) < in.S_MIN)
            In = I(1:in.S_MIN, n);
        end
        if (length(In) > in.S_MAX)
            In = I(1:in.S_MAX, n);
        end
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

        in.X(:, n) = 0;
        in.X(In, n) = coeff;
    end
    in.S_ESTIMATED = mean(in.SPARSITY);
else
    in.S_ESTIMATED = 1;
end
warning on;



end
