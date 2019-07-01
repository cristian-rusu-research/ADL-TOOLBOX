function [in] = itkrm_online(in)
in.IP = in.DICT'*in.Y;
absip = abs(in.IP);
signip = sign(in.IP);
[~,I] = sort(absip,1,'descend');
in.G = in.DICT'*in.DICT;

dnew = zeros(in.D,in.K);
avenapp = 0;
for n=1:in.N
    In = I(1:in.S,n);
    coeff = in.G(In,In)\in.IP(In,n);
    %%% if the line above creates a warning (instable)
    %%% recalculate stably but more expensively via pinv
    [~, msgid] = lastwarn;
    if ~isempty(msgid)
        if strcmp(msgid, 'MATLAB:nearlySingularMatrix') || ...
            strcmp(msgid, 'MATLAB:singularMatrix')
            coeff = pinv(gram(In,In))*ip(In,n);
        end
    end
    
    app = in.DICT(:,In)*coeff;
    res = in.Y(:,n) - app;
    enapp = app'*app;

    dnew(:,In) = dnew(:,In) + bsxfun(@times, repmat(res, 1, in.S), signip(In,n)');
    dnew(:,In) = dnew(:,In) + bsxfun(@times, in.DICT(:,In), absip(In,n)');
%     dnew(:,In) = dnew(:,In) + res*signip(In,n)';
%     dnew(:,In) = dnew(:,In) + in.DICT(:,In)*diag(absip(In,n));
    
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
