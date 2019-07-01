function [in] = plotdico(in)
if (isfield(in, 'FIG_DICT'))
    figure(in.FIG_DICT); imagesc(showdico(in.DICT)); axis off;
else
    in.FIG_DICT = figure; imagesc(showdico(in.DICT)); axis off;
end

if (isfield(in, 'FIG_DICT_TEXT'))
    title(['dictionary of ' in.FIG_DICT_TEXT ' , iter = ' num2str(in.ITER_CURRENT) ' , s = ' num2str(in.S) ' , K = ' num2str(in.K)]);
end

if (isfield(in, 'COHERENCE_MAX'))
    c = in.statistics(in.ITER_CURRENT).coherence;
    title(['iter = ' num2str(in.ITER_CURRENT) ' , s = ' num2str(in.S) ' , K = ' num2str(in.K) ' , coherence = '  num2str(c(end))]);
end

pause(0.5);

end
