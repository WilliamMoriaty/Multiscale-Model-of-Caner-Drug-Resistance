fig1=figure(1);
clf();
set(gcf,'Position',[488,531.4,675.4,230.6])

load('data.mat')
PC9_D=PC9_DRPs/100;
PC9_n=PC9_naive/100;

x=cumsum(PC9_n(:,2));
y=cumsum(PC9_D(:,2));

%% Using Naive data
phi0=0.063;phi1=0.25;
n=6;alpha=4.0;
eta=25;

subplot(1,2,1)
x=[PC9_n(:,1) x];

u=TurePar(alpha,phi0,phi1,n,eta);
y1=linspace(0,0.5,101);
[u1,y1]=hist(u,y1);
u1=cumsum(u1)/sum(u1);
plot(x(:,1),x(:,2),'Marker','pentagram',...
    'MarkerSize',12,'LineStyle','none','MarkerEdgeColor','none','MarkerFaceColor',[0.93,0.69,0.13])
hold on
plot(y1,u1,'LineWidth',1.2,'Color','k')

%title(strcat('\alpha=',num2str(alpha(1),4),',','\beta=',num2str(eta-alpha,4)))
title("Drug-naive cells (Naive)",'FontName','Helvetica','FontSize',13,'FontWeight','bold')
xlabel('CTP(x)')
ylabel('CDF')
set(gca,'FontName','Helvetica','FontSize',13,'FontWeight','bold','linewidth',1.2)

lgd=legend('data','fit');
lgd.FontWeight = 'bold';
% lgd.Position=[0.686 0.718 0.212 0.190];
lgd.Box='off';
lgd.ItemTokenSize = [10,6];

box off
%% change parameters for DRPs
alpha=6.0;

subplot(1,2,2)
y=[PC9_D(:,1) y];

u=TurePar(alpha,phi0,phi1,n,eta);
y2=linspace(0,1,201);
[u2,y2]=hist(u,y2);
u2=cumsum(u2)/sum(u2);
plot(y(:,1),y(:,2),'Marker','pentagram',...
    'MarkerSize',12,'LineStyle','none','MarkerEdgeColor','none','MarkerFaceColor',[0.93,0.69,0.13])
hold on
plot(y2,u2,'LineWidth',1.2,'Color','k')

ylim([0,1.2])
title("Drug-released persisters (DRPS)",'FontName','Helvetica','FontSize',13,'FontWeight','bold')
xlabel('CTP(x)')
ylabel('CDF')
set(gca,'FontName','Helvetica','FontSize',13,...
    'FontWeight','bold','linewidth',1.2,'ytick',0:0.2:1.2)

box off

%% iter algorithm for the distribution of epigentic state

function u=TurePar(alpha,phi0,phi1,n,eta)
U=zeros(101,2000);
U(1,:)=rand(1,2000);
for i=1:100
for k=1:2000
phi=phi0+phi1*(alpha*U(i,k))^n./(1+(alpha*U(i,k))^n);
a=eta*phi;
b=eta*(1-phi);
U(i+1,k)=betarnd(a,b);
end
end
u=U(end,:);
end


