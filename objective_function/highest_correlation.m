function [val] = highest_correlation(in, val)
if (isfield(in, 'SYNTHETIC_DICT'))
    val.highest_correlation = max(abs(in.DICT'*in.SYNTHETIC_DICT));
end
end
