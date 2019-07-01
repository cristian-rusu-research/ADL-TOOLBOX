function [in] = generate_repcand_ksvd_batch(in)

% [in.ITER_CURRENT,in.L]

%%% number of candidate iterations per dictionary iteration 
rep_iterations = ceil(log(in.D));
%%% corresponding number of training signals
repbatch_size = floor(size(in.Y,2)/rep_iterations);

%%%% corresponding threshold depending on risk (to do...)
in.REPLACEMENT_CANDIDATE_THRESHOLD = 1/sqrt(in.D);

% initialize random candidate atoms
repcand = randn(in.D, in.L);
repcand = bsxfun(@rdivide, repcand, sqrt(sum(repcand.^2)));
repfreq = zeros(in.L, 1);

allres = in.Y - in.DICT*in.X;

% train candidates on the residuals
for it = 1:rep_iterations
    res = allres(:,(it-1)*repbatch_size+1:it*repbatch_size);
    enres = sum(res.*res);
    absrepip = abs(repcand'*res);
    if in.L > 1
        [maxrepips, ind]=max(absrepip);
    else 
        ind = ones(1,repbatch_size);
        maxrepips = absrepip;
    end
    for j0=1:in.L
        % signals that use atom j0
        J = find(ind == j0);
        
        if (length(J) < 1)
    %         x = in.Y(:, randsample(in.N, 1));
    %         x = x/norm(x);
    %         in.DICT(:, j0) = x;
            continue;
        end
        
        [u, ~, ~] = svds(res(:, J), 1);
        repcand(:, j0)= u;
       
        if it == rep_iterations
            repfreq(j0) = sum(maxrepips(J).^2 >= enres(J)*in.REPLACEMENT_CANDIDATE_THRESHOLD);
        end
    end
 
    scale = sum(repcand.*repcand);
    repcand = bsxfun(@rdivide, repcand, sqrt(scale));
end

[repfreq,repind] = sort(repfreq, 'descend');
repcand = repcand(:, repind);

in.CANDIDATES = repcand;
in.CANDIDATE_SCORES = repfreq;

end
