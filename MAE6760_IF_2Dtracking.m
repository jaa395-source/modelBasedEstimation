%% Information Filter Example:
%   tracking of a pedestrian walking (e.g. phone based)
%   four states: x,y position and velocity
%   can select a different number of sensors
%
%   Note: uses plot_estimator.m, calculate_ellipse.m
%
% MAE 6760 Model Based Estimation
% Cornell University
% M Campbell
%
close all;clear all;
MAE6760startup; %can adjust font size, figure size at the bottom of this script
global MCcolors; %define colors as global with access 
rng(2);
%
%% Define 2D model and simulate
% system matrices: continuous 
Fc=[0 1 0 0;
    0 0 0 0;
    0 0 0 1;
    0 0 0 0]; %state: x,xdot,y,ydot
Gc=[0 0;1 0;0 0;0 1]; %process noise on both states x,y
Cc=[1 0 0 0;0 0 1 0]; %measure 2D position
nx=length(Fc);[~,nw]=size(Gc);[nz,~]=size(Cc);
%convert to discrete time
dt=0.1;
[sys]=c2d(ss(Fc,Gc,Cc,zeros(nz,nw)),dt,'zoh');
F=sys.A;G=sys.B;H=sys.C;

%Create truth via a person walking: 
%nominal case is a person walking straight, stopping, turning right
nt=100;
xv=[[0.1:0.1:1] ones(1,40) [1:-0.1:0.1] zeros(1,10) zeros(1,nt-70)]*2;
yv=[zeros(1,70) [0.1:0.1:1] ones(1,nt-70-10)]*2;
t=[0:dt:dt*(nt-1)];
%simulate to get position
xp=0;yp=0;
for k=1:(nt-1),
    xp(1,k+1) = xp(1,k)+dt*xv(k); %no process noise
    yp(1,k+1) = yp(1,k)+dt*yv(k); %no process noise
end
x_true=[xp;xv;yp;yv];

%Create 2D position GPS-like measurements 
R=eye(nz)*1^2;
v=sqrtm(R)*randn(nz,nt);
z=H*x_true + v;

%plot truth: x,y position over time
figure;
pp=plot(t,xp,'-',t,yp,':','Color',MCcolors.blue);
grid;
xlabel('time (sec)');ylabel('position');
legend('x position','y position','location','northwest');
axis([0 10 -1 11]);

%plot truth: birds-eye-view of the pedestrian position
figure;
pt=plot(xp,yp,'.','color',MCcolors.mag,'LineWidth',3); %true x,y position
hold on;
pp=plot(xp(1),yp(1),'>',xp(end),yp(end),'o','MarkerSize',10,'color',MCcolors.mag);
hold off;
axis([-1 11 -1 6]);
grid;
xlabel('x position (m)');ylabel('y position (m)');
legend('2D true position','start','stop','Location','Northwest');

%Initial conditions and process noise for all filters
x0=[0;0;0;0];
P0=diag([1^2 0.1^2 1^2 0.1^2]);
Q=eye(nw)*5; %non-zero Q a bit larger to cover unknown control


%% KF: Kalman Filter estimator
% pre-allocate memory
xhatu=zeros(nx,nt);xhatu(:,1)=x0; %define vector; initialize
Pu=zeros(nx,nx,nt);Pu(:,:,1)=P0; %define vector; initialize
xhatp=xhatu; %define vector; initialize
Pp=Pu; %define vector; initialize
for k=1:(nt-1),
    %predict
    xhatp(:,k+1) = F*xhatu(:,k);
    Pp(:,:,k+1) = F*Pu(:,:,k)*F' + G*Q*G';
    %Kalman Gain
    K = Pp(:,:,k+1)*H'*inv(H*Pp(:,:,k+1)*H' + R);
    %update
    xhatu(:,k+1) = xhatp(:,k+1) + K *(z(:,k+1) - H*xhatp(:,k+1));
    Pu(:,:,k+1) = (eye(nx)-K*H)*Pp(:,:,k+1)*(eye(nx)-K*H)' + K*R*K';
end

%plot KF: x,y position estimates over time
plot_estimator2(t,xhatu,Pu,x_true,z,[1 3]);
sgtitle('Kalman Filter (1 msmt)','FontWeight','bold','Fontsize',16);

%plot KF: birds-eye-view of the pedestrian position estimates
plot_birdseye(t,xhatu,Pu,x_true,[1 3]);
title('Kalman Filter (1 msmt)','FontWeight','bold');


%% IF: Information Filter estimator (1 msmt)
%IF initialization
Y0=inv(P0);
y0=Y0*x0;
yhatu=zeros(nx,nt);yhatu(:,1)=y0; %define vector; initialize
Yu=zeros(nx,nx,nt);Yu(:,:,1)=Y0; %define vector; initialize
%define physical state for plotting
xinfo=zeros(nx,nt);xinfo(:,1)=x0; %define vector; initialize
Pinfo=zeros(nx,nx,nt);Pinfo(:,:,1)=P0; %define vector; initialize
%inversion of system matrices for faster iteration
Fi=inv(F);Qi=inv(Q);Ri=inv(R);
for k=1:(nt-1),
    %predict
    Mi=Fi'*Yu(:,:,k)*Fi;
    Om=Mi*G*inv(Qi+G'*Mi*G);
    Yp(:,:,k+1) = (eye(nx)-Om*G')*Mi;   
    yhatp(:,k+1) = [(eye(nx)-Om*G')*Fi']*yhatu(:,k);
    %update
    yhatu(:,k+1) = yhatp(:,k+1) + H'*Ri*z(:,k+1);
    Yu(:,:,k+1) = Yp(:,:,k+1) + H'*Ri*H;
    %
    %pull out state and covariance for plotting
    Pinfo(:,:,k+1)=inv(Yu(:,:,k+1));
    xinfo(:,k+1) = Pinfo(:,:,k+1)*yhatu(:,k+1);
end

%plot IF (1 msmt): x,y position estimates over time
plot_estimator2(t,xinfo,Pinfo,x_true,z,[1 3]);
sgtitle('Information Filter (1 msmt)','FontWeight','bold','Fontsize',16);

%plot IF (1 msmt): birds-eye-view of the pedestrian position estimates
plot_birdseye(t,xinfo,Pinfo,x_true,[1 3]);
title('Information Filter (1 msmt)','FontWeight','bold');

%% ns msmt case
%Create additional measurements
ns=10; %total number of measurements
v_ns{1}=v;
z_ns{1}=z;
for j=2:ns,
    v_ns{j}=sqrtm(R)*randn(nz,nt);
    z_ns{j}=H*x_true + v_ns{j};
end

%% IF: Information Filter estimator (ns msmts)
%IF initialization
Y0=inv(P0);
y0=Y0*x0;
yhatu=zeros(nx,nt);yhatu(:,1)=y0; %define vector; initialize
Yu=zeros(nx,nx,nt);Yu(:,:,1)=Y0; %define vector; initialize
%define physical state for plotting
xinfo_ns=zeros(nx,nt);xinfo_ns(:,1)=x0; %define vector; initialize
Pinfo_ns=zeros(nx,nx,nt);Pinfo_ns(:,:,1)=P0; %define vector; initialize
%inversion of system matrices for faster iteration
Fi=inv(F);Qi=inv(Q);Ri=inv(R);
for k=1:(nt-1),
    %predict step
    Mi=Fi'*Yu(:,:,k)*Fi;
    Om=Mi*G*inv(Qi+G'*Mi*G);
    Yp(:,:,k+1) = (eye(nx)-Om*G')*Mi;   
    yhatp(:,k+1) = [(eye(nx)-Om*G')*Fi']*yhatu(:,k);
    %update step
    %find total info from all sensors first
    itot=zeros(nx,1);for j=1:ns,itot=itot+[H'*Ri*z_ns{j}(:,k+1)];end
    Itot=zeros(nx,nx);for j=1:ns,Itot=Itot+[H'*Ri*H];end
    %IF update
    yhatu(:,k+1) = yhatp(:,k+1) + itot;
    Yu(:,:,k+1) = Yp(:,:,k+1) + Itot;
    %
    %pull out state and covariance for plotting
    Pinfo_ns(:,:,k+1)=inv(Yu(:,:,k+1));
    xinfo_ns(:,k+1) = Pinfo_ns(:,:,k+1)*yhatu(:,k+1);
end

%plot IF (10 msmts): x,y position estimates over time
plot_estimator2(t,xinfo_ns,Pinfo_ns,x_true,z,[1 3]);
sgtitle('Information Filter (10 msmts)','FontWeight','bold','Fontsize',16);

%plot IF (10 msmts): birds-eye-view of the pedestrian position estimates
plot_birdseye(t,xinfo_ns,Pinfo_ns,x_true,[1 3]);
title('Information Filter (10 msmts)','FontWeight','bold');


%% --------- plotting functions
function plot_estimator2(t,xhat,Phat,x_true,z,ii);
%
figure('Position',[100 100 1600 600]);
tiledlayout(1,2,'TileSpacing','compact','Padding','tight');
nexttile;
plot_estimator(t,xhat(ii(1),:),Phat(ii(1),ii(1),:),x_true(ii(1),:),'error',z(1,:))
ylabel('x position (m)');
axis([0 10 -2 2]);
nexttile;
plot_estimator(t,xhat(ii(2),:),Phat(ii(2),ii(2),:),x_true(ii(2),:),'error',z(2,:))
ylabel('y position (m)');
axis([0 10 -2 2]);
%
end

function plot_birdseye(t,xhat,Phat,x_true,ii);
%
global MCcolors; %define colors as global with access 
figure;
pt=plot(x_true(ii(1),:),x_true(ii(2),:),'.','color',MCcolors.mag,'LineWidth',2); %true x,y position
hold on;
ph=plot(xhat(ii(1),:),xhat(ii(2),:),'-','color',MCcolors.blue); %state estimate
%plot 2-sigma uncertainty ellipsoids at a few positions along trajectory
iell=[2 10:10:length(t)];
for i=1:length(iell),
    [Xe,Ye] = calculate_ellipse(xhat(ii,iell(i)),Phat(ii,ii,iell(i)),2);
    pu(i)=plot(xhat(ii(1),iell(i)),xhat(ii(2),iell(i)),'x','color',MCcolors.blue,'LineWidth',1);
    ne=length(Xe);
    x1=[Xe((ne/2+1):ne);Xe(1:ne/2)];
    y1=[Ye((ne/2+1):ne);Ye(1:ne/2)];    
    pp(i)=patch(x1,y1,'b','facecolor',MCcolors.blue,'edgecolor',MCcolors.blue,'linewidth',1,'FaceAlpha',.1,'EdgeAlpha',1);

end
xlabel('x position (m)');ylabel('y position (m)');
hold off;
legend([ph pp(1) pt],'position estimate','error ellipse','truth','Location','Northwest');
%
end

%% --------- start-up items
function MAE6760startup(font_size);
%
% define colors for plotting
global MCcolors; %define colors as global with access 
MCcolors.red=[200,0,0]/255;
MCcolors.blue=[4,51,255]/255;
MCcolors.purple=[147,23,255]/255;
MCcolors.green=[0,160,0]/255;
MCcolors.orange=[253,128,8]/255;
MCcolors.mag=[255,64,255]/255;
MCcolors.cyan=[0,230,255]/255;
%
% define standard figure positioning and size
set(groot,'DefaultFigureUnits','pixels');
set(groot,'DefaultFigurePosition',[100 100 800 600]);
set(groot,'DefaultFigureWindowStyle','normal');  % Important
set(groot,'DefaultAxesFontSize',16);
set(groot,'DefaultAxesFontWeight','bold');
set(groot,'DefaultLineLineWidth',2);
%
end