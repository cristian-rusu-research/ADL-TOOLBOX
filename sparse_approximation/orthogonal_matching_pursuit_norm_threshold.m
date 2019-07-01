function [in] = orthogonal_matching_pursuit_norm_threshold(in)

if (~strcmp(in.SPARSITY_MODE, 'ADAPTIVE'))
    error('This implementation of matching pursuit sparse recovery needs adaptive sparsity!');
end

if (~isfield(in, 'SPARSITY_APPROXIMATION_ERROR'))
    in.SPARSITY_APPROXIMATION_ERROR = 1/5;
    warning(['Do not use this function if you do not have a reasonable way to estimate the noise level! Default SPARSITY_APPROXIMATION_ERROR is set to ' num2str(in.SPARSITY_APPROXIMATION_ERROR)]);
end

for n = 1:in.N 
    
    y = in.Y(:,n);
    normy = norm(y);
    supp = [];
    coeff = [];
    GramInvNew = zeros(in.S_MAX);

    res = y;
    normres = normy;
    k = 0;
    
    absip = abs(in.DICT'*res);
    [maxip, i] = max(absip);
    
    while ((k < in.S_MAX) && (normres > in.SPARSITY_APPROXIMATION_ERROR*normy))  || (k < in.S_MIN)
        i = i(1);
        newatom = in.DICT(:, i);
        supp = [supp i];

        if k == 0,
            %Initializing everything
            GramInvNew(1,1) = 1;
            k = 1;
            GramInv = 1/norm(newatom)^2;
        else
            %Update the Inverse of the Gram Matrix
            Q = GramInv * dico_supp' * newatom;
            m = 1 / (newatom' * newatom - ...
                 newatom' * dico_supp * GramInv * dico_supp' * newatom);
            GramInvNew(1:k, 1:k) = GramInv + Q * Q' * m;
            GramInvNew(1:k, k+1) =  - Q * m;
            GramInvNew(k+1, 1:k) =  - Q' * m;
            GramInvNew(k+1, k+1) = m;
            k = k + 1;
            GramInv = GramInvNew(1:k , 1:k);
        end	
        
        dico_supp = in.DICT(:,supp);
        
        coeff = GramInv * (dico_supp' * y);
        res = y - dico_supp * coeff;
        normres = norm(res);
        
        absip = abs(in.DICT'*res);
        [maxip,i] = max(absip);
    
    end
    
    in.X(:, n) = 0;
    in.X(supp, n) = coeff;
end

end
