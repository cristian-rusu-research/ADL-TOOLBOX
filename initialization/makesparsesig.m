function sigN=makesparsesig(dico,N,S,T,b,rho,bfix)

% syntax: sigN = makesparsesig(dico,N,S,T,b,rho)
%
% makes sparse signals according to the 
% model described in Table 1 of
% 'Local Identification of Overcomplete Dictionaries'
% arXiv:1401.6354
%
% input:
% dico... generating dictionary
% N... number of sparse signals to create
% S... effective sparsity - number of strong coefficients - default 1
% T... number of non-zero coefficients T>=S - default T=S
% b... decay parameter for coefficients - default b=0
% rho... noiselevel - default - rho=0
%
% output:
% sigN... d x N matrix with N sparse signals as columns
%
% Karin Schnass 24.01.14

if (nargin<2)	
    disp('syntax: sigN=makesparsesig(dico,N,S,T,b,rho)');
    sigN=[];		
    return;
end

inputflaw=false;

if (nargin<3)
    S=1;
end
if round(S)<1
    S=max([round(abs(S)),1]);
    inputflaw=true;
end
if (nargin<4)
    T=S;
end

if T<S
    inputflaw=true;
    T=S;
end

if (nargin<5)
    b=0;
end	

if nargin < 6
   rho=0;
end

if (nargin<7)
    bfix=1;
end

if (b<0) || (b>1)
    inputflaw=true;
    b=0;
end

[d, K]=size(dico);

if S>d-1
   inputflaw=true;
   S=d-1;
end

if T>K;
   inputflaw=true;
   T=K;
end


if inputflaw==true
   disp('warning, strange input parameters, used: [d,S,T,b,rho]=');
   [d,S,T,b,rho]
end

sigN=[];

for n=1:N
    if bfix == 1
        beta=1-b;
        
    else
        beta=1-b *rand(1,1);
    end
    x1toS=sqrt(1./S)*beta.^[1:S]';
    x1toSsign=2*round(rand(S,1))-1;
    x1toS=x1toS.*x1toSsign;
    %norm(x1toS)
    if T > S
        xSp1toT=randn(T-S,1);
        xSp1toT=xSp1toT* sqrt(1-norm(x1toS)^2)/norm(xSp1toT);
        x1toT=[x1toS; xSp1toT];
    else
        x1toT=x1toS/norm(x1toS);
    end
    p=randperm(K);
    sig=dico(:,p(1:T))*x1toT;
    if (rho > 0)
        noise = rho*randn(d,1);
        sig=(sig+noise)/sqrt(1+noise'*noise);
    end     

    sigN=[sigN,sig];  
    
end

