function [in] = default_prune_max_coherent_atoms(in)

if (in.K == 1)
    return;
end

%%%% add code to not go below K_min-L

if (~isfield(in, 'COHERENT_ATOM_REPLACEMENT_MODE'))
	in.COHERENT_ATOM_REPLACEMENT_MODE = 'm';
	warning('COHERENT_ATOM_REPLACEMENT_MODE is set to ''m''');
end

% keep track of atoms to prune
prune = [];

% compute Gram matrix
holgram = abs(in.DICT'*in.DICT - eye(in.K));
maxcorr = max(max(holgram));
while (maxcorr > in.COHERENCE_MAX) && (length(prune) < in.K-in.K_MIN + in.L)
    [rowi, coli] = find(holgram == maxcorr, 1, 'first');
    h = sign(in.DICT(:, rowi)'*in.DICT(:, coli));
    
    if ismember(in.COHERENT_ATOM_REPLACEMENT_MODE,['c', 'd'])
        % compare most recent score
        if in.ATOM_SCORES(end, rowi) > in.ATOM_SCORES(end, coli)
            newatom = in.DICT(:, rowi);
        else
            newatom = in.DICT(:, coli);
        end
    elseif in.COHERENT_ATOM_REPLACEMENT_MODE == 'm'
        % correctly signed sum weighted according to most recent score
        newatom = in.ATOM_SCORES(end, rowi)*in.DICT(:, rowi) + h*in.ATOM_SCORES(end, coli)*in.DICT(:, coli);
    else 
        % correctly signed sum 
        newatom = in.DICT(:, rowi)+ h*in.DICT(:, coli);
    end
    
    if in.ATOM_SCORES(end, rowi) < in.ATOM_SCORES(end, coli)
        dummy=coli;
        coli = rowi;
        rowi = dummy;
    end
    
    in.DICT(:, rowi) = newatom;
    if norm(newatom) > 0
        in.DICT(:, rowi)= newatom/norm(newatom);
    end
    in.ATOM_SCORES(end, rowi) = in.ATOM_SCORES(end, rowi) + in.ATOM_SCORES(end, coli);
    %in.X(rowi,:) = in.X(rowi,:) + h*in.X(coli,:);
    prune = [prune, coli];
    
    holgram(:, [rowi, coli]) = 0;
    holgram([rowi, coli], :) = 0;
    maxcorr = max(max(holgram));
end

in.DICT(:, prune) = [];
in.K = size(in.DICT, 2);
in.X(prune, :) = [];
in.ATOM_SCORES(:, prune) = [];

end
