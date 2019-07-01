function y = checkMembers(in, list)

for i = 1:length(list)
    if (~isfield(in, list(i)))
        y = 0; return;
    end
end

y = 1;
