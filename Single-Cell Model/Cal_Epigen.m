function S = Cal_Epigen(alpha0,alpha1,eta,c,N)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

n=100;
x=linspace(1/n,1,n);
z=linspace(1/n,1,n);
P=zeros(n);kappa=0.063;
zeta0=0.25;
for i=1:n
for j=1:n
P(i,j)=inheritance(x(j),z(i),eta,kappa,alpha0,zeta0);
end
end

P1=zeros(n);
for i=1:n
P1(i,:)=P(i,:)./sum(P(i,:));
end

P=zeros(n);
kappa=0.063;
zeta1=0.25;
for i=1:n
for j=1:n
P(i,j)=inheritance(x(j),z(i),eta,kappa,alpha1,zeta1);
end
end

P2=zeros(n);
for i=1:n
P2(i,:)=P(i,:)./sum(P(i,:));
end


y=rand(1,n);
y=y./sum(y);
% S=zeros(1500,100);
S=zeros(N,100);
S(1,:)=y;dtau=4e-2;
beta0=0.4;beta1=1.0;
a1=8;a2=9;a3=15;
beta=beta0+beta1.*(a1*x+(a2*x).^6)./(1+(a3*x).^6);
for i=2:N
if i<301 || c(i)==0
S_in=2*(beta.*S(i-1,:)*P1);
S_out=S(i-1,:).*(beta+sum(beta.*S(i-1,:)));
S(i,:)=S(i-1,:)+dtau*(S_in-S_out);
else 
S_in=2*(beta.*S(i-1,:)*P2);
S_out=S(i-1,:).*(beta+sum(beta.*S(i-1,:)));
S(i,:)=S(i-1,:)+dtau*(S_in-S_out);
end
end

end