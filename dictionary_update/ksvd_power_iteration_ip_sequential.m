function [in] = ksvd_power_iteration_ip_sequential(in)
absip = abs(in.IP);

dnew = zeros(in.D,in.K);
avenapp = 0;
for n=1:in.N
    In = in.SUPPORTS(:, n);
    
    app = in.DICT(:,In)*in.X(In, n);
    res = in.Y(:,n) - app;
    enapp = app'*app;
    
    dnew(:,In) = dnew(:,In) + bsxfun(@times, repmat(res, 1, in.S), in.IP(In, n)');
    dnew(:,In) = dnew(:,In) + bsxfun(@times, in.DICT(:,In), (absip(In, n).*in.IP(In, n))');
%     dnew(:,In)=dnew(:,In) + real(res*ip(In,n)');
%     dnew(:,In)=dnew(:,In) + dico(:,In)*diag(absip(In,n).*ip(In,n));
    
    avenapp = avenapp + enapp;
end

avenapp = avenapp/in.N;

scale = sum(dnew.*dnew);
nonzero = find(scale > 0.001*avenapp/in.D);
% iszero = find(scale <= 0.001*avenapp/in.D);
iszero = setdiff(1:in.K, nonzero);

% dnew(:,nonzero) = dnew(:,nonzero)*diag(1./sqrt(scale(nonzero)));
dnew(:,nonzero) = bsxfun(@rdivide, dnew(:,nonzero), sqrt(scale(nonzero)));
in.DICT(:,nonzero) = dnew(:,nonzero);

end
