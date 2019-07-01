function [in] = ksvd_batch(in)
%% update of K-SVD (update both DICT and X)

p = randperm(in.K);
for j = 1:in.K
    j0 = p(j);
        
    J = find(in.X(j0,:));
    
    if (length(J)<1)
%         x = in.Y(:, randsample(in.N, 1));
%         x = x/norm(x);
%         in.DICT(:, j0) = x;
        continue;
    end
    
    support = [1:j0-1 j0+1:in.K];
    ERj0 = in.Y(:, J) - in.DICT(:, support)*in.X(support, J);

%     options.tol = 1e-2;
    [u, s, v] = svds(ERj0, 1, 'L');
    in.DICT(:, j0) = u;
    in.X(j0, J) = s*v';
end

end
