function y = norms(x, dim)
y = sqrt( sum( x .* conj( x ), dim ) );
