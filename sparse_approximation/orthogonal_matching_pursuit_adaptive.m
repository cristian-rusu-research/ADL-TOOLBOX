function [in] = orthogonal_matching_pursuit_adaptive(in)

if (~strcmp(in.SPARSITY_MODE, 'ADAPTIVE'))
    error('This implementation of matching pursuit sparse recovery needs adaptive sparsity!');
end

if (~isfield(in, 'SPARSITY_NOISE_RISK'))
    in.SPARSITY_NOISE_RISK = 1/4;
    warning(['SPARSITY_NOISE_RISK is set to ' num2str(in.SPARSITY_NOISE_RISK)]);
end

threshold = sqrt(2*(log(2*in.K) - log(in.SPARSITY_NOISE_RISK))/in.D);

for n = 1:in.N 
    
    y = in.Y(:, n);
    normy = norm(y);
    supp = [];
    coeff = [];
    GramInvNew = zeros(in.S_MAX);

    res = y;
    normres = normy;
    k = 0;
    
    absip = abs(in.DICT'*res);
    [maxip, i] = max(absip);
    
    Smax = min(in.K,in.S_MAX);
    
    while ((maxip > threshold*normres) && (k < Smax)) || (k < in.S_MIN),
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
