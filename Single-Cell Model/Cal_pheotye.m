function h = Cal_pheotye(t,h1,par,c,S)
% Calculate phenotype
% x: number
% h: m*n matrix,
% m:time step,n:phenotype vector

%% parameter
n=1001;m=100;
h=zeros(m,n);

deltat=par.deltat;

x=linspace(1/m,1,m);
y=linspace(0,1.0,n);
%y=zeros(1,n);



if t==0
%% initial state
for i=1:m    
h(i,:)=exp(-(y-0.2).^2/par.Delta);
h(i,:)=h(i,:)/(par.dy*sum(h(i,:)));
end
else

%% explicity iteration
for i=1:m
ii=ceil(t/100);
if (t<30001)
v=par.v*(par.y0(1)-y);
gamma2=0;
gamma=gamma2*c*(y<par.yth);
elseif c>0
gamma2=1.0;
v=par.v*c*x(i);
v=v*(par.y0(2)-y);
gamma=gamma2*c*(y<par.yth);
else
v=par.v*x(i)*(par.y0(1)-y);
gamma2=0;
gamma=gamma2*c*(y<par.yth);
end
h(i,2:n-1)=h1(i,1:n-2).*(deltat*v(1:n-2)/(2*par.dy)+par.Delta*deltat/par.dy^2)...
    +(1-2*par.Delta*deltat/par.dy^2).*h1(i,2:n-1)...
    +h1(i,3:n).*(-deltat*v(3:n)/(2*par.dy)+par.Delta*deltat/par.dy^2)...
    -deltat*h1(i,2:n-1).*(gamma(2:n-1)-sum(sum(gamma.*h1.*S(ii,:)')));
%% boundary
h(i,1)=0;
h(i,n)=0;
%% normal one
h(i,:)=h(i,:)/(sum(h(i,:)));
end

end

end