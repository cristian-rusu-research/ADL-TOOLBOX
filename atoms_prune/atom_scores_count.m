function [in] = atom_scores_count(in)

in.ATOM_SCORES = [in.ATOM_SCORES; sum(abs(sign(in.X)), 2)'];
    
if (size(in.ATOM_SCORES, 1) > in.ATOM_SCORES_MEMORY)
    in.ATOM_SCORES = in.ATOM_SCORES(2:end, :);
end

if (~isfield(in, 'COHERENCE_MAX'))
    plot(in.ATOM_SCORES(end, :)); title(['iteration ' num2str(in.ITER_CURRENT)]);
    xlabel('atom index'); ylabel('importance score');
    drawnow
    pause(0.5);
end

% maxscore = max(in.ATOM_SCORES);
% if (maxscore > 1)
%     in.THRESHOLD_REMOVE_UNUSED_ATOMS = [in.THRESHOLD_REMOVE_UNUSED_ATOMS min(1, maxscore/2)];
% else
%     % safeguard for ill-chosen number of observations
%     % (compared to number of training signals)
%     % and data that is not sparse
%     in.THRESHOLD_REMOVE_UNUSED_ATOMS = [in.THRESHOLD_REMOVE_UNUSED_ATOMS maxscore/sqrt(in.D)];
% end

end
