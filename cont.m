close all
clear
clc

R1=0.1; 
R2=1;
th=[0:1:180] * pi/180;
x1=R1*sin(th);
y1=R1*cos(th);
x2=R2*sin(th);
y2=R2*cos(th);
x3=[0 0];
y3=[R1 R2];
x4=[0 0];
y4=[-R1 -R2];

KN=77; LN=77; MN=77;

for i=1:100
 map(i,1)=0.0;
 map(i,2)=0.0;
 map(i,3)=0.0;
end
%%%%%%%%%%%% 

load fort.11
x=reshape(fort(:,1),201,181);
y=reshape(fort(:,2),201,181);

load fort.12
g=reshape(fort(:,1),201,181);
h=reshape(fort(:,2),201,181);

load fort.13
e=reshape(fort(:,1),201,181);
f=reshape(fort(:,2),201,181);

%%%%%%%%%%%%

emin=min(min(e));
emax=max(max(e));

fmin=min(min(f));
fmax=max(max(f));

gmin=min(min(g));
gmax=max(max(g));

hmin=min(min(h));
hmax=max(max(h));

%%%%%% e
figure
axis equal
axis off
hold
plot(x1,y1,'k',x2,y2,'k',x3,y3,'k',x4,y4,'k')
contour(x,y,e,[0:emax/10:emax],'r')
contour(x,y,e,[emin:-emin/10:0],'b')
%set(gcf, 'Renderer', 'zbuffer')
colormap(map)
print -dpdfcrop angular.pdf

%%%%%% f
figure
axis equal
axis off
hold
plot(x1,y1,'k',x2,y2,'k',x3,y3,'k',x4,y4,'k')
contour(x,y,f,[0:fmax/10:fmax],'r')
contour(x,y,f,[fmin:-fmin/10:0],'b')
%set(gcf, 'Renderer', 'zbuffer')
colormap(map)
print -dpdfcrop circulation.pdf

%%%%%% g
figure
axis equal
axis off
hold
plot(x1,y1,'k',x2,y2,'k',x3,y3,'k',x4,y4,'k')
contour(x,y,g,[0:gmax/10:gmax],'r')
contour(x,y,g,[gmin:-gmin/10:0],'b')
%set(gcf, 'Renderer', 'zbuffer')
colormap(map)
print -dpdfcrop toroidal.pdf

%%%%%% h
figure
axis equal
axis off
hold
plot(x1,y1,'k',x2,y2,'k',x3,y3,'k',x4,y4,'k')
contour(x,y,h,[0:hmax/10:hmax],'r')
contour(x,y,h,[hmin:-hmin/10:0],'b')
%set(gcf, 'Renderer', 'zbuffer')
colormap(map)
print -dpdfcrop poloidal.pdf

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

th=[0:1:360]*pi/180;
x5=R1*sin(th);
y5=R1*cos(th);
x6=R2*sin(th);
y6=R2*cos(th);

load fort.16
XE=reshape(fort(:,1),2*MN,KN);
YE=reshape(fort(:,2),2*MN,KN);

load fort.17
B1=reshape(fort(:,1),2*MN,KN);
B2=reshape(fort(:,2),2*MN,KN);
B3=reshape(fort(:,3),2*MN,KN);

load fort.18
U1=reshape(fort(:,1),2*MN,KN);
U2=reshape(fort(:,2),2*MN,KN);
U3=reshape(fort(:,3),2*MN,KN);

XE=[XE(end,:);XE];
YE=[YE(end,:);YE];
B1=[B1(end,:);B1];
B2=[B2(end,:);B2];
B3=[B3(end,:);B3];
U1=[U1(end,:);U1];
U2=[U2(end,:);U2];
U3=[U3(end,:);U3];

B1max=max(max(B1));
B1min=min(min(B1));
B2max=max(max(B2));
B2min=min(min(B2));
B3max=max(max(B3));
B3min=min(min(B3));
U1max=max(max(U1));
U1min=min(min(U1));
U2max=max(max(U2));
U2min=min(min(U2));
U3max=max(max(U3));
U3min=min(min(U3));

figure
axis equal
axis off
hold
plot(x5,y5,'k',x6,y6,'k')
contour(XE,YE,U1,[0:U1max/10:U1max],'r')
contour(XE,YE,U1,[U1min:-U1min/10:0],'b')
%set(gcf, 'Renderer', 'zbuffer')
colormap(map)

figure
axis equal
axis off
hold
plot(x5,y5,'k',x6,y6,'k')
contour(XE,YE,U2,[0:U2max/10:U2max],'r')
contour(XE,YE,U2,[U2min:-U2min/10:0],'b')
%set(gcf, 'Renderer', 'zbuffer')
colormap(map)

figure
axis equal
axis off
hold
plot(x5,y5,'k',x6,y6,'k')
contour(XE,YE,U3,[0:U3max/10:U3max],'r')
contour(XE,YE,U3,[U3min:-U3min/10:0],'b')
%set(gcf, 'Renderer', 'zbuffer')
colormap(map)
print -dpdfcrop uphi.pdf

figure
axis equal
axis off
hold
plot(x5,y5,'k',x6,y6,'k')
contour(XE,YE,B1,[0:B1max/10:B1max],'r')
contour(XE,YE,B1,[B1min:-B1min/10:0],'b')
%set(gcf, 'Renderer', 'zbuffer')
colormap(map)

figure
axis equal
axis off
hold
plot(x5,y5,'k',x6,y6,'k')
contour(XE,YE,B2,[0:B2max/10:B2max],'r')
contour(XE,YE,B2,[B2min:-B2min/10:0],'b')
%set(gcf, 'Renderer', 'zbuffer')
colormap(map)

figure
axis equal
axis off
hold
plot(x5,y5,'k',x6,y6,'k')
contour(XE,YE,B3,[0:B3max/10:B3max],'r')
contour(XE,YE,B3,[B3min:-B3min/10:0],'b')
%set(gcf, 'Renderer', 'zbuffer')
colormap(map)
print -dpdfcrop bphi.pdf

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load fort.19
xm=reshape(fort(:,1),LN,KN);
ym=reshape(fort(:,2),LN,KN);
helicity=reshape(fort(:,3),LN,KN);
helicitymin=min(min(helicity));
helicitymax=max(max(helicity));
figure
axis equal
axis off
hold
plot(x1,y1,'k',x2,y2,'k',x3,y3,'k',x4,y4,'k')
contour(xm,ym,helicity,[0:helicitymax/10:helicitymax],'r')
contour(xm,ym,helicity,[helicitymin:-helicitymin/10:0],'b')
%set(gcf, 'Renderer', 'zbuffer')
colormap(map)
