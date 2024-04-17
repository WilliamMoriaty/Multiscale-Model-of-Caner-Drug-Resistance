nn=1;

% dd=8:19;
cc=linspace(1/20,1,20);
for ii=1:nn

%parameter
z=1;
T=150000;
N=round(T/100); % N=T/100;
hh=zeros(100,1001);
H=zeros(100,1001,N);
F=zeros(100,1001,N);

alpha0=4.0;alpha1=alpha0+(1/0.3)*0.9;
eta=25;
%% calculate epigenetic CTP
c=zeros(7201,1);
c(300:end)=cc(ii);
%dd=8:19;
% 25 cycle = 1 day
% d=dd(ii);
% for i=d:d:288
% c(300+(i-d)*25:301+(i)*25,1)=0.9;
% c(300+7*25+(i-d)*25:301+(i)*25,1)=0;
% end



% end
%load('Low_dose.mat')

S = Cal_Epigen(alpha0,alpha1,eta,c,N);
par=struct;
par.v=1.5; %2.5
par.Delta=1e-3;
par.dy=1/1000;
par.y0=[0.2 0.8];
par.yth=0.5;
par.deltat=0.4*par.dy^2/par.Delta;

for t=0:T
h = Cal_pheotye(t,hh,par,c(z),S);
hh=h;
if mod(t,100)==0 && t>0
H(:,:,z)=h;
z=z+1;
end

end


for i=1:N
for j=1:100
F(j,1:1001,i)=S(i,j)*H(j,1:end,i);
end
end

% FF=F(:,:,end);
%path=strcat("Drug_alpha",num2str(ii),"_eta",num2str(jj),".mat");
%path='Drug_Low_dose_HH_eta.mat';
path=strcat("Drug_Dose","_",num2str(ii),".mat");
save(path,'S','H','F');

end



