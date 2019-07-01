function [in] = atom_scores_weighted_count(in)

if (~isfield(in, 'TOTAL_DATASET_ENERGY'))
    in.TOTAL_DATASET_ENERGY = norm(in.Y,'fro')^2;
end

in.ATOM_SCORES = [in.ATOM_SCORES; sum(abs(in.X).^2,2)'/in.TOTAL_DATASET_ENERGY];
    
if (size(in.ATOM_SCORES, 1) > in.ATOM_SCORES_MEMORY)
    in.ATOM_SCORES = in.ATOM_SCORES(2:end, :);
end

if (isfield(in, 'FIG_IMPORTANCE'))
    figure(in.FIG_IMPORTANCE);
else
    in.FIG_IMPORTANCE = figure;
end
    plot(in.ATOM_SCORES(end, :)); title(['iteration ' num2str(in.ITER_CURRENT)]);
    xlabel('atom index'); ylabel('importance score');
    drawnow

if (isfield(in, 'FIG_DICT'))
    figure(in.FIG_DICT);
else
    in.FIG_DICT = figure;
end

     imagesc(showdico(in.DICT)); axis off;
     title(['dictionary atoms, average sparsity is ' num2str(length(find(in.X))./in.N)]);

end
