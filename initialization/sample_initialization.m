function [in] = sample_initialization(in)
%% update only the dictionary randomly
in.DICT = in.Y(:, randsample(1:in.N, in.K));
in.DICT = bsxfun(@rdivide, in.DICT, sqrt(sum(in.DICT.^2)));

%% set the sparse representation to zero, they are updated first in the iterative process
% in.X = sparse(in.K, in.N);
in.X = zeros(in.K, in.N);
in.SUPPORTS = zeros(in.S, in.N);

end
