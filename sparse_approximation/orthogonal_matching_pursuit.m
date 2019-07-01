function [in] = orthogonal_matching_pursuit(in)

if (~strcmp(in.SPARSITY_MODE, 'FIXED'))
    error('This implementation of orthogonal matching pursuit needs fixed sparsity!');
end

in.X = omp(in.DICT'*in.Y, in.DICT'*in.DICT, in.S);

for n = 1:in.N
    in.SUPPORTS(:, n) = find(in.X(:,n));
end

end
