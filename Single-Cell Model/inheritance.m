function p=inheritance(x,z,eta,kappa,alpha,zeta)
% p(x,z):probability of inheritant function
% x:daughter cells;y: parent cells
phi=kappa+zeta*(alpha*z)^6./(1+(alpha*z)^6);
a=eta*phi;
b=eta*(1-phi);

p=x^(a-1)*(1-x)^(b-1)/beta(a,b);


end