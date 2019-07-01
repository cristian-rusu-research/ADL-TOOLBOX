function y = rule_remove_coherent_atoms(in)

if (in.ITER_CURRENT >= round(log(in.D)))
    y = true;
else
    y = false;
end

end
