function [in] = add_atoms(in)

if (in.K == in.K_MAX)
    return;
end

newatomscores = in.THRESHOLD_REMOVE_UNUSED_ATOMS*ones(size(in.ATOM_SCORES, 1), 1);

for j0 = 1:in.L
    if (in.K == in.K_MAX)
        break;
    end
        
    if (in.CANDIDATE_SCORES(j0) > 1 || in.K < in.K_MIN )
        if (isfield(in, 'COHERENCE_MAX'))
            rcoh = max(abs(in.DICT'*in.CANDIDATES(:, j0)));
            if (rcoh < in.COHERENCE_MAX || in.K < in.K_MIN)
                in.DICT = [in.DICT in.CANDIDATES(:,j0)];
                in.K = in.K + 1;
                in.X = [in.X; zeros(1,size(in.Y,2))];
                in.ATOM_SCORES = [in.ATOM_SCORES, newatomscores];
            end
        else
            in.DICT = [in.DICT in.CANDIDATES(:,j0)];
            in.K = in.K + 1;
            in.X = [in.X; zeros(1,size(in.Y,2))];
            in.ATOM_SCORES = [in.ATOM_SCORES, newatomscores];
        end
    end
end

end
