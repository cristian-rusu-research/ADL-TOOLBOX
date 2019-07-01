function y = rule_remove_unused_atoms(in)

if (in.ITER_CURRENT >= 2*in.ATOM_SCORES_MEMORY && mod(in.ITER_CURRENT, 1) == 0 && in.ITER_CURRENT <= in.ITER_MAXIMUM - 10) && in.ITER_CURRENT < (in.ITER_MAXIMUM - in.ATOM_SCORES_MEMORY)
    y = true;
else
    y = false;
end
    
end
