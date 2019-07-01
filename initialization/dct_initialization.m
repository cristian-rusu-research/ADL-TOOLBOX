function [in] = dct_initialization(in)
%% set dictionary to DCT samples
if (in.K <= in.D)
    in.DICT = dctmtx(in.D);
    in.DICT = in.DICT(:, 1:in.K);
else
    if (~isfield(in, 'DCT_LINEAR_COMBINATIONS'))
        warning('DCT_LINEAR_COMBINATIONS set to 3.');
        in.DCT_LINEAR_COMBINATIONS = 3;
    end
    
    C = sparse(in.D, in.K);
    for i = 1:in.K
        support = randsample(1:in.D, in.DCT_LINEAR_COMBINATIONS);
        C(support, i) = randn(in.DCT_LINEAR_COMBINATIONS, 1);
    end
    in.DICT = dct(C);
end

in.DICT = bsxfun(@rdivide, in.DICT, sqrt(sum(in.DICT.^2)));

%% set the sparse representation to zero, they are updated first in the iterative process
% in.X = sparse(in.K, in.N);
in.X = zeros(in.K, in.N);
in.SUPPORTS = zeros(in.S, in.N);

end
