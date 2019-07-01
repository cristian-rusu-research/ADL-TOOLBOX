function [in] = orthogonal(in)
%% Orthogonal dictionary learning

if (in.K ~= in.D)
    error('Orthogonal dictionaries can be created only square dictionaries, K = D!');
end

[U, ~, V] = svd(in.Y*in.X');
in.DICT = U*V';
