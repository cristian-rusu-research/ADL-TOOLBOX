function [val] = frobenius_norm_squared_dictionary_size_coherence(in, val)
if (~isfield(in, 'TOTAL_DATASET_ENERGY'))
    in.TOTAL_DATASET_ENERGY = norm(in.Y, 'fro')^2;
end

val.frobenius_norm_squared = norm(in.Y - in.DICT*in.X, 'fro')^2/in.TOTAL_DATASET_ENERGY;
val.coherence = max(max(abs(in.DICT'*in.DICT - eye(in.K))));
val.dictionary_size = in.K;
end
