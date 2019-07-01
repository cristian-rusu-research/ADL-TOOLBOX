function [in] = method_of_directions(in)
%% MOD update

in.DICT = (in.Y*in.X')/(in.X*in.X');
in.DICT = bsxfun(@rdivide, in.DICT, sqrt(sum(in.DICT.^2)));
