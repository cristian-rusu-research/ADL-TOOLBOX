function [in] = atom_scores_count_above_threshold(in)
app = in.DICT*in.X;
res = in.Y - app;

enres = sum(res.^2)*in.TAU;
enapp = sum(app.^2);
threshold = (enres + enapp)/in.D;

in.ATOM_SCORES = [in.ATOM_SCORES; sum(abs(in.X).^2 > threshold, 2)'];
    
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
