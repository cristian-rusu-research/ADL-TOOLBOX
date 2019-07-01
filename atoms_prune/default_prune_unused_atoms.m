function [in] = default_prune_unused_atoms(in)

if (in.K == 1)
    return;
end

%%%% add code to not go below K_min-L

maxprune = min(min(round(in.D/5),floor(in.K/2)), in.K - in.K_MIN + in.L);

notused = find(max(in.ATOM_SCORES, [], 1) < in.THRESHOLD_REMOVE_UNUSED_ATOMS);

if length(notused) > maxprune
    [~, I] = sort(max(in.ATOM_SCORES, [], 1), 'ascend');
    notused = I(1:maxprune);
end

in.DICT(:, notused) = [];
in.K = size(in.DICT, 2);
in.X(notused, :) = [];
in.ATOM_SCORES(:, notused) = [];

end
