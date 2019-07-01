function [val] = frobenius_norm(in, val)
val.frobenius_norm_squared = norm(in.Y - in.DICT*in.X, 'fro');
end
