%% MAE 6760 Model Based Estimation
% Cornell University
% M Campbell
%
% Homework #3
% Problem #2: Kalman Filter
%   Two aircraft tracking estimator in 2D as they turn
%   Plotting uses plot_openloop.m and plot_estimator.m
%
clear all; close all
MAE6760startup; %can adjust font size, figure size at the bottom of this script
global MCcolors; %define colors as global with access 
rng(101);
 
%% Define system model (CT and DT) for aircraft A 
%continuous time model for aircraft A
Omega_A = 0.045; %turn velocity for A in rad/s
%state vector: [East pos, East vel, North pos, North vel]
nx=4;
A_A=[0 1 0 0;0 0 0 -Omega_A;0 0 0 1;0 Omega_A 0 0]; 
ix=1;iy=3; %indices of the 2D position
%input forces (East, North)
B_A=[0 0;1 0;0 0;0 1]; 
nw=2;
%discrete time model 
dt = 0.5; %sec
[F_A,G_A]=c2d(A_A,B_A,dt);
% sensor output models: 2D position states
H_A = [1 0 0 0; 
       0 0 1 0];
nz=2;
% build discrete time model for aircraft A
DTsys_A = ss(F_A,G_A,H_A,zeros(nz,nw),dt); 
%model directional wind disturbances
Q_A = 10*[2.0 0.05;    
         0.05 0.5];

%% simulate aircraft A with process noise and plot 
tvec = 0:dt:100;nk=length(tvec);
x0_A = [0,85*cos(pi/4),0,-85*sin(pi/4)]';
w_A = sqrtm(Q_A)*randn(nw,nk);
[~,~,x_A] = lsim(DTsys_A,w_A,tvec,x0_A);x_A=x_A'; %make into row vector
z_A_clean = H_A*x_A;
%
%simulate measurements for aircraft A from a tracking station
R_A = [20 0.05; 
      0.05 20];
v_A = sqrtm(R_A)*randn(nz,nk);
z_A = H_A*x_A+v_A;

%plot trajectory for aircraft A for some context
figs(1)=figure('Position',[100 100 800 600]);
plot_openloop_2Daircraft(x_A);
title('2D trajectory for aircraft A')


%% Part a) for aircraft A: generate measurements and create a Kalman Filter
%goal: KF to track/localize aircraft A
%
%use model DTsys_A for this problem
%use measurements z_A for this problem
%

F=DTsys_A.A;G=DTsys_A.B;H=H_A;
nx=length(F);
Q=Q_A;R=R_A;z=z_A;

P0_val = 1;
X0_val = 1;
P0 = diag(P0_val*ones(nx,1));

xhatp = zeros(nx,nk);xhatp(:,1) = x0_A;
xhatu = zeros(nx,nk);xhatu(:,1) = x0_A;
Pp = zeros(nx,nx,nk);Pp(:,:,1) = P0;
Pu = zeros(nx,nx,nk);Pu(:,:,1) = P0;
n = nk;

for k=1:(n-1)
    % Preidct
    xhatp(:,k+1) = F*xhatu(:,k);
    Pp(:,:,k+1) = F*Pu(:,:,k)*F' + G*Q*G';

    % Kalman Gain
    K = Pp(:,:,k+1)*H'*inv(H*Pp(:,:,k+1)*H' + R);

    % Update
    xhatu(:,k+1) = xhatp(:,k+1) +K*(z(:,k+1) - H*xhatp(:,k+1));
    Pu(:,:,k+1) = (eye(nx) - K*H)*Pp(:,:,k+1)*(eye(nx) - K*H)' + K*R*K';
end


%Generate plot of East and North position error estimates and 2-sigma bounds
%assumes t is (1 x nk), xhatu is (nx x nk), Pu is (nx x nx x nk), xtrue is (nx x nk) 
%
figs(2)=figure;
ti(2)=tiledlayout(1,2,'TileSpacing','compact','Padding','tight');
nexttile;
plot_estimator(tvec,xhatu(ix,:),Pu(ix,ix,:),x_A(ix,:),'error',z_A(1,:));
ylabel('North error estimate {\it{e_N(t)}}');
nexttile;
plot_estimator(tvec,xhatu(iy,:),Pu(iy,iy,:),x_A(iy,:),'error',z_A(2,:));
ylabel('East error estimate {\it{e_E(t)}}');
title(ti(2),'(a): Baseline Kalman Filter for aircraft A','fontweight','bold');


%% Part b): same as Part a), but design 95% measurement gate
%goal: KF to track/localize aircraft A with measurement gate
%
%use model DTsys_A for this problem
%use measurements z_A for this problem
%
%variables to generate for display/plots:
% Nrej: number of measurement rejections (subset of total nk)
% Trej: times [tvec(k+1)] when msmt is rejected 
% Erej: error [z(:,k+1)-H*x_A(:,k+1)] when msmt is rejected

%msmt gating parameters
Pgate = 0.95;alpha = 1-Pgate;
Lam0  = chi2inv(Pgate, nz);
Nrej=0;Irej=[];

F=DTsys_A.A;G=DTsys_A.B;H=H_A;
nx=length(F);
Q=Q_A;R=R_A;z=z_A;

%%YOUR CODE HERE

%uncomment lines below to output the percent of msmts rejected
% disp('(b) percent of msmts rejected:');
% disp(Nrej/nk*100)

%Generate plot of East and North position error estimates and 2-sigma bounds
%assumes t is (1 x nk), xhatu is (nx x nk), Pu is (nx x nx x nk), xtrue is (nx x nk) 
figs(3)=figure;
ti(3)=tiledlayout(1,2,'TileSpacing','compact','Padding','tight');
nexttile;
plot_estimator(tvec,xhatu(ix,:),Pu(ix,ix,:),x_A(ix,:),'error',z_A(1,:));
ylabel('North error estimate {\it{e_N(t)}}');
hold on
%uncomment line below to add msmt rejections to the plot
%plot(Trej,Erej(1,:),'mo','DisplayName','Msmt Rejection');
hold off;
nexttile;
plot_estimator(tvec,xhatu(iy,:),Pu(iy,iy,:),x_A(iy,:),'error',z_A(2,:));
ylabel('East error estimate {\it{e_E(t)}}');
hold on
%uncomment line below to add msmt rejections to the plot
%plot(Trej,Erej(2,:),'mo','DisplayName','Msmt Rejection');
hold off;
title(ti(3),'(b): Kalman Filter with 95% msmt gating for aircraft A','fontweight','bold');


%% Part c): Estimate Process Noise Covariance Q to make KF consistent
%goal: use hypothesis tests to tune process noise covariance Q 
%
%use model DTsys_A for this problem
%load simulated measurements from tile (these are for aircraft A, but labeled z_c for Part c
%use measurements z_c for this problem
%use states x_c for the plots on problem
%
%variables to generate for display/plots:
% Lam: innovations test statistic lambda (1 x nk)
% LamF: Kalman Filter test statistic lambda (1 x nk)
% Nrej: number of measurement rejections (subset of total nk)
% Trej: times [tvec(k+1)] when msmt is rejected 
% Erej: error [z(:,k+1)-H*x_c(:,k+1)] when msmt is rejected
% NFrej: number of times the filter is inconsistent (subset of total nk)
% TFrej: times [tvec(k+1)] when the filter is inconsistent 
% EFrej: error [z(:,k+1)-H*x_c(:,k+1)] when the filter is inconsistent
%
%load measurement vector z_c and state x_c 
% for aircraft A, simulated with different Q
load ACdata 
%

F=DTsys_A.A;G=DTsys_A.B;H=H_A;
nx=length(F);
Q=Q_A;R=R_A;z=z_c;

%msmt gating parameters
Pgate = 0.95;alpha = 1-Pgate;
Lam0  = chi2inv(Pgate, nz);

%filter consistency parameters
PF = 0.95;alpha = 1-PF;
win=10;%time window
Blow=chi2inv(alpha/2,10*nz)/win;  %low filter threshold
Bhigh=chi2inv(1-alpha/2,10*nz)/win;  %high filter threshold
Lam=zeros(nk,1);
LamF=zeros(nk,1);

%%YOUR CODE HERE

%
%uncomment lines below to output the percent of msmts rejected
% disp('(c) percent of msmts rejected:');
% disp(Nrej/nk*100)
%uncomment lines below to output the percent of time filter is inconsistent
% disp('(c) percent of time filter is inconsistent:');
% disp(NFrej/nk*100)
%
%Generate plot of East and North position error estimates and 2-sigma bounds
%assumes t is (1 x nk), xhatu is (nx x nk), Pu is (nx x nx x nk), xtrue is (nx x nk) 
figs(4)=figure;
ti(4)=tiledlayout(1,2,'TileSpacing','compact','Padding','tight');
nexttile;
plot_estimator(tvec,xhatu(ix,:),Pu(ix,ix,:),x_c(ix,:),'error',z_c(1,:));
ylabel('North error estimate {\it{e_N(t)}}');
hold on
%uncomment the line below to add msmt rejections to the plot
plot(Trej,Erej(1,:),'mo','DisplayName','Msmt Rejection');
hold off;
nexttile;
plot_estimator(tvec,xhatu(iy,:),Pu(iy,iy,:),x_c(iy,:),'error',z_c(2,:));
ylabel('East error estimate {\it{e_E(t)}}');
hold on
%uncomment the line below to add msmt rejections to the plot
plot(Trej,Erej(2,:),'mo','DisplayName','Msmt Rejection');
hold off;
title(ti(4),'(c): Kalman Filter with Msmt Gate and Filter Integrity for aircraft A','fontweight','bold');

figs(5)=figure('Position',[100 100 800 600]);
plot(tvec,Lam,':','Color',MCcolors.mag);
hold on;
plot(tvec((win+1):end),LamF((win+1):end),'-','Color',MCcolors.purple);
plot(tvec,Blow*ones(1,nk),'--','Color',MCcolors.blue);
plot(tvec,Bhigh*ones(1,nk),'--','Color',MCcolors.blue);
hold off;
ylabel('filter integrity')
xlabel('time \it{t} (sec)')
title('(c) Kalman Filter integrity test statistic')
legend('inn test statistic \lambda','KF test statistic \lambda_{10}^{KF}','lower bound','upper bound')


%% Part d): One joint KF for aircraft A and aircraft B
%goal: estimate X=[x_A;x_B] given Z=[z_A;z_B]
%
%continuous time model for aircraft B
Omega_B = -0.045; %turn velocity for V in rad/s
%state vector: [East pos, East vel, North pos, North vel]
A_B=[0 1 0 0;0 0 0 -Omega_B;0 0 0 1;0 Omega_B 0 0]; 
%input forces (East, North)
B_B=[0 0;1 0;0 0;0 1]; 
%discretize 
[F_B,G_B]=c2d(A_B,B_B,dt);
% sensor output models: 2D position states
H_B = [1 0 0 0; 
       0 0 1 0];
%model directional wind disturbances
Q_B = 10*[2.0 0.05;    
         0.05 0.5];
%
% build discrete time model for AC A
DTsys_B = ss(F_B,G_B,H_B,zeros(nz,nw),dt); 

%simulate aircraft B with process noise 
x0_B = [4000,85*cos(pi/4),3200,-85*sin(pi/4)]';
w_B = sqrtm(Q_B)*randn(nw,nk);
[~,~,x_B] = lsim(DTsys_B,w_B,tvec,x0_B);x_B=x_B'; %make into row vector
z_B_clean = H_B * x_B;

%simulate measurements for aircraft A from a tracking station
R_B = [20 0.05; 
      0.05 20];
v_B = sqrtm(R_B)*randn(nz,nk);
z_B = z_B_clean+v_B;
%

%plot trajectories for both aircraft A and B for some context
figs(6)=figure('Position',[100 100 800 600]);
plot_openloop_2Daircraft(x_A,x_B);
title('2D trajectory for aircraft A and B')
%
%use both models DTsys_A,DTsys_B for this problem
%use measurements z_A,z_B for this problem
%
F=blkdiag(DTsys_A.A,DTsys_B.A);
G=blkdiag(DTsys_A.B,DTsys_B.B);
H=blkdiag(H_A,H_B);
nX=length(F);
Q=blkdiag(Q_A,Q_B);
R=blkdiag(R_A,R_B);
Z=[z_A;z_B];

%%YOUR CODE HERE

%Generate plot of East and North position error estimates and 2-sigma bounds
%assumes t is (1 x nk), xhatu is (nx x nk), Pu is (nx x nx x nk), xtrue is (nx x nk) 
figs(7)=figure;
ti(7)=tiledlayout(1,2,'TileSpacing','compact','Padding','tight');
nexttile;
plot_estimator(tvec,xhatu(ix,:),Pu(ix,ix,:),x_A(ix,:),'error',z_A(1,:));
ylabel('North error estimate {\it{e_N(t)}}');
nexttile;
plot_estimator(tvec,xhatu(iy,:),Pu(iy,iy,:),x_A(iy,:),'error',z_A(2,:));
ylabel('East error estimate {\it{e_E(t)}}');
title(ti(7),'(d): Joint Kalman Filter: aircraft A errors','fontweight','bold');
figs(8)=figure;
ti(8)=tiledlayout(1,2,'TileSpacing','compact','Padding','tight');
nexttile;
plot_estimator(tvec,xhatu(nx+ix,:),Pu(nx+ix,nx+ix,:),x_B(ix,:),'error',z_B(1,:));
ylabel('North error estimate {\it{e_N(t)}}');
nexttile;
plot_estimator(tvec,xhatu(nx+iy,:),Pu(nx+iy,nx+iy,:),x_B(iy,:),'error',z_B(2,:));
ylabel('East error estimate {\it{e_E(t)}}');
title(ti(8),'(d): Joint Kalman Filter: aircraft B errors','fontweight','bold');


%% Part e) same as Part d), but uses a 2D relative range sensor
%goal: estimate X=[X_A;X_B] given Z=[z_A;z_B]
%      --> repeat (d), but for a 2D range measurement
%
%use both models DTsys_A,DTsys_B for this problem
%use measurements Zr for this problem
%
Rr = [10 0.15; 
      0.15 10];
Hr = [1 0 0 0 -1 0  0 0;
         0 0 1 0  0 0 -1 0]; 
vr = sqrtm(Rr)*randn(2,nk);
Zr = [x_A(1,:)-x_B(1,:);x_A(3,:)-x_B(3,:)]+vr;

%use both DTsys_A and DTsys_B for this problem
F=blkdiag(DTsys_A.A,DTsys_B.A);
G=blkdiag(DTsys_A.B,DTsys_B.B);
H=Hr;
nX=length(F);
Q=blkdiag(Q_A,Q_B);
R=Rr;
Z=Zr;

%%YOUR CODE HERE

%Generate plot of East and North position error estimates and 2-sigma bounds
%assumes t is (1 x nk), xhatu is (nx x nk), Pu is (nx x nx x nk), xtrue is (nx x nk) 
figs(9)=figure;
ti(9)=tiledlayout(1,2,'TileSpacing','compact','Padding','tight');
nexttile;
plot_estimator(tvec,xhatu(ix,:),Pu(ix,ix,:),x_A(ix,:),'error',z_A(1,:));
ylabel('North error estimate {\it{e_N(t)}}');
nexttile;
plot_estimator(tvec,xhatu(iy,:),Pu(iy,iy,:),x_A(iy,:),'error',z_A(2,:));
ylabel('East error estimate {\it{e_E(t)}}');
title(ti(10),'(e): Joint KF with Relative Range: aircraft A errors','fontweight','bold');
figs(10)=figure;
ti(10)=tiledlayout(1,2,'TileSpacing','compact','Padding','tight');
nexttile;
plot_estimator(tvec,xhatu(nx+ix,:),Pu(nx+ix,nx+ix,:),x_B(ix,:),'error',z_B(1,:));
ylabel('North error estimate {\it{e_N(t)}}');
nexttile;
plot_estimator(tvec,xhatu(nx+iy,:),Pu(nx+iy,nx+iy,:),x_B(iy,:),'error',z_B(2,:));
ylabel('East error estimate {\it{e_E(t)}}');
title(ti(10),'(e): Joint KF with Relative Range: aircraft B errors','fontweight','bold');


%% Plotting functions
function plot_openloop_2Daircraft(x_A,x_B);
% x_A is the state of aircraft A [East pos, East vel, North pos, North vel]
% optional: x_B is the state of aircraft B [East pos, East vel, North pos, North vel]
%
global MCcolors;
%plot aircraft A
p1=plot(x_A(1,:),x_A(3,:),'-','color',MCcolors.blue);
hold on;
plot(x_A(1,1),x_A(3,1),'b>','markersize',10);
plot(x_A(1,end),x_A(3,end),'pentagram','Color','b','markersize',10);
hold off
ylabel('North');xlabel('East');
grid;
%
if nargin<2, %only aircraft A  
    legend([p1],'aircraft A');
else, %overlay aircraft B
    hold on;
    p2=plot(x_B(1,:),x_B(3,:),':','color',MCcolors.red);
    hold on;
    plot(x_B(1,1),x_B(3,1),'>','color',MCcolors.red,'markersize',10);
    plot(x_B(1,end),x_B(3,end),'pentagram','Color',MCcolors.red,'markersize',10);
    hold off
    legend([p1 p2],'aircraft A','aircraft B');
end
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
set(groot,'DefaultFigurePosition',[100 100 1600 600]);
set(groot,'DefaultFigureWindowStyle','normal');  % Important
set(groot,'DefaultAxesFontSize',16);
set(groot,'DefaultAxesFontWeight','bold');
set(groot,'DefaultLineLineWidth',2);
%
end