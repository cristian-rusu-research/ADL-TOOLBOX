function in = remove_mean(in)
%% remove signal's means and return
in.Y = bsxfun(@minus, in.Y, mean(in.Y));
in.TOTAL_DATASET_ENERGY = norm(in.Y, 'fro')^2;
