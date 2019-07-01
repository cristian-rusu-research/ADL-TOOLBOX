function y = rule_add_atoms(in)

m = round(log(in.D));
if (in.ITER_CURRENT >= m && in.ITER_CURRENT < (in.ITER_MAXIMUM - 5*m) && mod(in.ITER_CURRENT, 1) == 0)
    y = true;
else
    y = false;
end
    
end
