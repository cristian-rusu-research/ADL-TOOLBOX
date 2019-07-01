function [in] = ploterror(in)
if (isfield(in, 'FIG_ERROR'))
    figure(in.FIG_ERROR); plot([in.statistics(1:in.ITER_CURRENT).frobenius_norm_squared]);
else
    in.FIG_ERROR = figure; plot(in.statistics.frobenius_norm_squared);
end

if (isfield(in, 'FIG_ERROR_TEXT'))
    title(['dictionary of ' in.FIG_ERROR_TEXT ' , iter = ' num2str(in.ITER_CURRENT) ' , s = ' num2str(in.S) ' , K = ' num2str(in.K)]);
end

xlabel('iteration');
ylabel('representation error (frobenius norm squared)');
pause(0.5);

end
