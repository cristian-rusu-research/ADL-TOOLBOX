function [in] = thresholding(in)

if (~strcmp(in.SPARSITY_MODE, 'FIXED'))
    error('This implementation of thresholding sparse recovery needs fixed sparsity!');
end

ip = in.DICT'*in.Y;
absip = abs(ip);
[~, I] = sort(absip, 1, 'descend');
G = in.DICT'*in.DICT;

%[in.S in.K]

for n=1:in.N
    In = I(1:in.S,n);
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
    in.SUPPORTS(:, n) = In;
end

end
