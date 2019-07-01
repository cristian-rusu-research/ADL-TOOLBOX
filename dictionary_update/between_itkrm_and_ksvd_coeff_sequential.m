function [in] = between_itkrm_and_ksvd_coeff_sequential(in)
dnew = zeros(in.D,in.K);
% avenapp = 0;
for n=1:in.N
    In = in.SUPPORTS(:, n);
    
    app = in.DICT(:,In)*in.X(In, n);
    res = in.Y(:,n) - app;
    enapp = app'*app;
    
    weight = in.X(In, n);
    dnew(:,In) = dnew(:,In) + bsxfun(@times, repmat(res, 1, in.S), ((sign(weight).*abs(weight)).^(in.P-1))');
    dnew(:,In) = dnew(:,In) + bsxfun(@times, in.DICT(:,In), ((sign(weight).*abs(weight)).^in.P)');
%     weight = coeff;
%     dnew(:,In)=dnew(:,In) + real(res*(sign(weight).*abs(weight).^(p-1))');
%     dnew(:,In)=dnew(:,In) + dico(:,In)*diag(sign(weight).*abs(weight).^p);
    
%     avenapp = avenapp + enapp;
end

% avenapp = avenapp/in.N;

scale = sum(dnew.*dnew);
% nonzero = find(scale > 0.001*avenapp/in.D);
% iszero = find(scale <= 0.001*avenapp/in.D);
% iszero = setdiff(1:in.K, nonzero);

% dnew(:,nonzero) = dnew(:,nonzero)*diag(1./sqrt(scale(nonzero)));
% dnew(:,nonzero) = bsxfun(@rdivide, dnew(:,nonzero), sqrt(scale(nonzero)));
% in.DICT(:,nonzero) = dnew(:,nonzero);
in.DICT = bsxfun(@rdivide, dnew, sqrt(scale));

end
