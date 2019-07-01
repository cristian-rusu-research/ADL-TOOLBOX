function [val] = coherence(in, val)
val.coherence = max(max(abs(in.DICT'*in.DICT - eye(in.K))));
end
