function in = default_generate_synthetic_dataset(in)

if ~isfield(in, 'N')
    in.N = 20*in.SYNTHETIC_DATA_DIMENSION;
end

in.Y = zeros(in.SYNTHETIC_DATA_DIMENSION, in.N);
Sl = max(1, in.SYNTHETIC_DATA_SPARSITY - in.SYNTHETIC_DATA_SPARSITY_RANGE);
Ninito4 = ceil(in.N/4);

in.Y(:, 1:Ninito4) = makesparsesig(in.SYNTHETIC_DICT, Ninito4, Sl, Sl, in.SYNTHETIC_DATA_DECAY, in.SYNTHETIC_DATA_NOISE_LEVEL, 0);
in.Y(:, Ninito4+1:2*Ninito4) = makesparsesig(in.SYNTHETIC_DICT, Ninito4, in.SYNTHETIC_DATA_SPARSITY + in.SYNTHETIC_DATA_SPARSITY_RANGE, in.SYNTHETIC_DATA_SPARSITY + in.SYNTHETIC_DATA_SPARSITY_RANGE, in.SYNTHETIC_DATA_DECAY, in.SYNTHETIC_DATA_NOISE_LEVEL, 0);
in.Y(:, 2*Ninito4+1:in.N) = makesparsesig(in.SYNTHETIC_DICT, in.N - 2*Ninito4, in.SYNTHETIC_DATA_SPARSITY, in.SYNTHETIC_DATA_SPARSITY, in.SYNTHETIC_DATA_DECAY, in.SYNTHETIC_DATA_NOISE_LEVEL, 0);
p = randperm(in.N);
in.Y = in.Y(:,p);
if in.SYNTHETIC_DATA_OUTLIERS > 0
    p = randperm(in.N);
    Nout = round(in.N*in.SYNTHETIC_DATA_OUTLIERS);
    Yout = randn(in.SYNTHETIC_DATA_DIMENSION, Nout)/sqrt(in.SYNTHETIC_DATA_DIMENSION);
    in.Y(:, p(1:Nout)) = Yout;
end

if in.SYNTHETIC_DATA_NORMALIZE > 0
    in.Y = bsxfun(@rdivide, in.Y, sqrt(sum(in.Y.^2)));
end
