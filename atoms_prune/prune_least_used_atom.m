function [in] = prune_least_used_atom(in)

if (in.K == 1)
    return;
end

% average atom scores
[~, I] = min(in.ATOM_SCORES(end, :));

in.DICT(:, I(1)) = [];
in.K = in.K - 1;
in.X(I(1), :) = [];

% reset scores
in.ATOM_SCORES(:, I(1)) = [];

end
