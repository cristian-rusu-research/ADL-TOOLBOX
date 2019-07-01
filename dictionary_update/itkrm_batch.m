function [in] = itkrm_batch(in)
ip = in.DICT'*in.Y;
absip = abs(ip);
signip = sign(ip);
clear ip;

dnew = zeros(in.D, in.K);

app = in.DICT*in.X;
res = in.Y - app;
avenapp = norm(app, 'fro');

for j0=1:in.K
    % signal that use atom j0
    J = find(in.X(j0, :));
    
    if (length(J) < 1)
%         x = in.Y(:, randsample(in.N, 1));
%         x = x/norm(x);
%         in.DICT(:, j0) = x;
        continue;
    end
    
    dnew(:, j0) = res(:, J)*signip(j0, J)';
    dnew(:, j0) = dnew(:, j0) + in.DICT(:, j0)*sum(absip(j0, J));
%     dnew(:,In) = dnew(:,In) + res*signip(In,n)';
%     dnew(:,In) = dnew(:,In) + in.DICT(:,In)*diag(absip(In,n));
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
