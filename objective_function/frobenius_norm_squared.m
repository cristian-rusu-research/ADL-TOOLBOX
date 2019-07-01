function [val] = frobenius_norm_squared(in, val)
if (~isfield(in, 'TOTAL_DATASET_ENERGY'))
    in.TOTAL_DATASET_ENERGY = norm(in.Y, 'fro')^2;
end

val.frobenius_norm_squared = norm(in.Y - in.DICT*in.X, 'fro')^2/in.TOTAL_DATASET_ENERGY;
end
