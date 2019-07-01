function [in] = random_initialization(in)
%% update only the dictionary randomly
in.DICT = randn(in.D, in.K);
in.DICT = bsxfun(@rdivide, in.DICT, sqrt(sum(in.DICT.^2)));

%% set the sparse representation to zero, they are updated first in the iterative process
% in.X = sparse(in.K, in.N);
in.X = zeros(in.K, in.N);
end
