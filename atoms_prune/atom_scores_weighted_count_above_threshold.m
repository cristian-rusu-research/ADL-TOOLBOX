function [in] = atom_scores_weighted_count_above_threshold(in)
app = in.DICT*in.X;
res = in.Y - app;

totalenergy = sum(sum(in.Y.*in.Y));

enres = sum(res.^2)*in.TAU;
enapp = sum(app.^2)/in.D;
threshold = enres + enapp;

aux = abs(in.X).^2;

in.ATOM_SCORES = [in.ATOM_SCORES; sum(aux.*(aux > threshold), 2)'];
    
if (size(in.ATOM_SCORES, 1) > in.ATOM_SCORES_MEMORY)
    in.ATOM_SCORES = in.ATOM_SCORES(2:end, :);
end

if (~isfield(in, 'COHERENCE_MAX'))
    plot(in.ATOM_SCORES(end, :)); title(['iteration ' num2str(in.ITER_CURRENT)]);
    xlabel('atom index'); ylabel('importance score');
    drawnow
else
    plot(out.statistics.coherence); title(['iteration ' num2str(in.ITER_CURRENT)]);
    xlabel('iteration'); ylabel('coherence');
    drawnow
end


end
