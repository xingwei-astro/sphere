clear; close all; clc
lu=50; mu=20; lb=50; mb=20;

load ekin.dat
l=reshape(ekin(:,1),lu,mu+1);
m=reshape(ekin(:,2),lu,mu+1);
ek=reshape(ekin(:,3),lu,mu+1);
ek=log10(ek);
load emag.dat
l=reshape(emag(:,1),lb,mb+1);
m=reshape(emag(:,2),lb,mb+1);
em=reshape(emag(:,3),lb,mb+1);
em=log10(em);
figure
subplot(1,2,1)
contour(l,m,ek,'showtext','on')
%colorbar
xlabel('degree $l$','fontsize',15,'interpreter','latex')
ylabel('order $m$','fontsize',15,'interpreter','latex')
title('contours of $log_{10}EK$','fontsize',15,'interpreter','latex')
set(gca, 'OuterPosition', [0, 0, 0.5, 0.5])
subplot(1,2,2)
contour(l,m,em,'showtext','on')
%colorbar
xlabel('degree $l$','fontsize',15,'interpreter','latex')
ylabel('order $m$','fontsize',15,'interpreter','latex')
title('contours of $log_{10}EM$','fontsize',15,'interpreter','latex')
set(gca, 'OuterPosition', [0.5, 0, 0.5, 0.5])
print -dpsc spec1.ps

load ekin-l.dat
load ekin-m.dat
load emag-l.dat
load emag-m.dat
figure
semilogy(ekin_l(:,1),ekin_l(:,2),'-k',ekin_m(:,1),ekin_m(:,2),'-r',emag_l(:,1),emag_l(:,2),'-g',emag_m(:,1),emag_m(:,2),'-b')
legend('EK-l','EK-m','EM-l','EM-m','fontsize',15)
print -dpsc spec2.ps
