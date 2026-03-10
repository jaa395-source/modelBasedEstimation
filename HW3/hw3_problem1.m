%% MAE 6760 Model Based Estimation
% Cornell University
% M Campbell
%
% Homework #3
% Problem #1: Kalman Filter
%   Five mass-spring-damper system with initial condition and disturbance
%
clear all; close all
MAE6760startup; %can adjust font size, figure size at the bottom of this script

%% Define system model: continuous and then discrete
k1=1;k2=1;k3=1;k4=1;k5=1;
c1=0.01;c2=0.01;c3=0.01;c4=0.01;c5=0.01;
m1=1;m2=1;m3=1;m4=1;m5=1;
%System matrix, with states q1,q1dot,...
nx=10;
Fc=zeros(nx,nx);
Fc([2 4 6 8 10],[1 3 5 7 9])=...
 [-(k1+k2)/m1 (k2)/m1 0 0 0; 
   (k2)/m2 -(k2+k3)/m2 (k3)/m2 0 0;
   0 (k3)/m3 -(k3+k4)/m3 (k4)/m3 0;
   0 0 (k4)/m4 -(k4+k5)/m4 (k5)/m4;
   0 0 0 (k5)/m5 -(k5)/m5];      
Fc([2 4 6 8 10],[2 4 6 8 10])=...
 [-(c1+c2)/m1 (c2)/m1 0 0 0; 
   (c2)/m2 -(c2+k3)/m2 (c3)/m2 0 0;
   0 (c3)/m3 -(c3+k4)/m3 (c4)/m3 0;
   0 0 (c4)/m4 -(c4+c5)/m4 (c5)/m4;
   0 0 0 (c5)/m5 -(c5)/m5]; 
Fc([1 3 5 7 9],[2 4 6 8 10])=eye(5);
%disturbance w on the 5th mass
nw=1;
Gc=zeros(nx,1);Gc(nx,1)=1/m5;
nz=2;
Hc=[1 0 zeros(1,8); %position of the 1st mass
    zeros(1,8) 1 0];%position of the 5th mass
%CT model
CTsys=ss(Fc,Gc,Hc,zeros(nz,nw));
%DT model
dt=0.05;
DTsys = c2d(CTsys,dt);

%% Simulate the dynamic system
%initial condition on the position 5th mass
x0=zeros(nx,1);x0(9)=1; %displace the 5th mass
t=[0:dt:50];nk=length(t);
rng(101) %sets the same random number generator seed
Qsim=0.01;
w=randn(nk,1)*sqrt(Qsim);
[~,tsim,x_no_w] = lsim(DTsys,w*0,t,x0); %response to just the initial condition
x_no_w=x_no_w'; %make into row vector
[~,tsim,x_true] = lsim(DTsys,w,t,x0);   %response to the initial condition and disturbance
x_true=x_true'; %make into row vector

%Create noisy measurements
i1=1;
Rq1=0.001;vq1=sqrt(Rq1)*randn(1,nk);
z_q1=x_true(i1,:)+vq1;
i5=9;
Rq5=0.001;vq5=sqrt(Rq5)*randn(1,nk);
z_q5=x_true(i5,:)+vq5;

%plots of the states and measurements for the 1st/5th masses
figs(1)=figure;
ti(1)=tiledlayout(1,2,'TileSpacing','compact','Padding','tight');
nexttile;
i1=1;
plot_openloop(t,x_true(i1,:),z_q1,x_no_w(i1,:))
ylabel('position of the 1st mass');xlabel('time (sec)');
nexttile;
i2=9;
plot_openloop(t,x_true(i2,:),z_q5,x_no_w(i2,:))
ylabel('position of the 5th mass');xlabel('time (sec)');
title(ti(1),'simulated response: 1st and 5th masses','fontweight','bold');


%% Part (a): two Kalman Filters using z_q1, z_q5
%KF using z_q1
F=DTsys.A;G=DTsys.B;Q=Qsim;
H=DTsys.C(1,:);
R=Rq1;z=z_q1;

P0_vals = [];
for i = 1:nx
    P0_vals(i) = cov(x_no_w(i,:));
end

x0(1) = 0;
x0(9) = 1; 
P0 = diag(P0_vals);
%Generate plot of position error estimates for the 1st and 5th masses and 2-sigma bounds
%assumes t is (1 x nk), xhatu is (nx x nk), Pu is (nx x nx x nk), x_true is (nx x nk) 
%measurement of position of 5th mass

xhatp = zeros(nx,nk);xhatp(:,1) = x0;
xhatu = zeros(nx,nk);xhatu(:,1) = x0;
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


figs(2)=figure;
ti(2)=tiledlayout(1,2,'TileSpacing','compact','Padding','tight');
nexttile;
plot_estimator(t,xhatu(i1,:),Pu(i1,i1,:),x_true(i1,:),'error',z);
ylabel('error position q1');
ylim([-0.3 0.3]);
nexttile;
plot_estimator(t,xhatu(i2,:),Pu(i2,i2,:),x_true(i2,:),'error');
ylabel('error position q5');
ylim([-2 2]);
title(ti(2),'(a): measurement z_{q1}','fontweight','bold');

% KF using z_q5
F=DTsys.A;G=DTsys.B;Q=Qsim;
H=DTsys.C(2,:);
R=Rq5;z=z_q5;


%Generate plot of position error estimates for the 1st and 5th masses and 2-sigma bounds
%assumes t is (1 x nk), xhatu is (nx x nk), Pu is (nx x nx x nk), x_true is (nx x nk) 
%measurement of position of 5th mass

xhatp = zeros(nx,nk);xhatp(:,1)= x0;
xhatu = zeros(nx,nk);xhatu(:,1)= x0;
Pp = zeros(nx,nx,nk);Pp(:,:,1)= P0;
Pu = zeros(nx,nx,nk);Pu(:,:,1)= P0;
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
figs(3)=figure;
ti(3)=tiledlayout(1,2,'TileSpacing','compact','Padding','tight');
nexttile;
plot_estimator(t,xhatu(i1,:),Pu(i1,i1,:),x_true(i1,:),'error');
ylabel('error position q1');
ylim([-1 1]);
nexttile;
plot_estimator(t,xhatu(i2,:),Pu(i2,i2,:),x_true(i2,:),'error',z);
ylabel('error position q5');
ylim([-0.3 0.3]);
title(ti(3),'(a): measurement z_{q5}','fontweight','bold');


%% Part (b): predicted vs updated covariance using z_q5
%KF using z_q5
F=DTsys.A;G=DTsys.B;Q=Qsim;
H=DTsys.C(2,:);
R=Rq5;z=z_q5;

xhatp = zeros(nx,nk);xhatp(:,1)= x0;
xhatu = zeros(nx,nk);xhatu(:,1)= x0;
Pp = zeros(nx,nx,nk);Pp(:,:,1)= P0;
Pu = zeros(nx,nx,nk);Pu(:,:,1)= P0;
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

figs(4)=figure('Position',[100 100 800 600]);
semilogy(t,[squeeze(Pp(i5,i5,:))],'b:',t,[squeeze(Pu(i5,i5,:))],'r-');
axis([0 50 1E-5 1E-3]);
grid
legend('predicted','updated');
ylabel('estimate');xlabel('time (sec)');
title('(b) predicted vs updated');  


%% Part (c): Steady State error covariance
F=DTsys.A;G=DTsys.B;Q=Qsim;
H=DTsys.C(2,:);
R=Rq5;z=z_q5;

Pz_sim = (Pu(:,:,end));
[Pz_ss,K_ss,L_ss] = idare(F',H',G*Q*G',R);


%simulated measurement covariance from the KF at the final time:
disp('The simulated measurement covariance at the final time is:')
disp([Pz_sim(5,5)]);

%steady state measurement covariance:
disp('The steady state measurement covariance is:')
disp([Pz_ss(5,5)]);


%% Part (d): Kalman Filter initialized with steady state error covariance
%KF using z_q5
F=DTsys.A;G=DTsys.B;Q=Qsim;
H=DTsys.C(2,:);
R=Rq5;z=z_q5;


xhatp = zeros(nx,nk);xhatp(:,1)= x0;
xhatu = zeros(nx,nk);xhatu(:,1)= x0;
Pp = zeros(nx,nx,nk);Pp(:,:,1)= Pz_ss;
Pu = zeros(nx,nx,nk);Pu(:,:,1)= Pz_ss;
n = nk;
for k=1:(n-1)
    % Preidct
    xhatp(:,k+1) = F*xhatu(:,k);
    Pp(:,:,k+1) = F*Pu(:,:,k)*F' + G*Q*G';

    % Kalman Gain
    K = K_ss';

    % Update
    xhatu(:,k+1) = xhatp(:,k+1) +K*(z(:,k+1) - H*xhatp(:,k+1));
    Pu(:,:,k+1) = (eye(nx) - K*H)*Pp(:,:,k+1)*(eye(nx) - K*H)' + K*R*K';
end

%Generate plot of position error estimates for the 1st and 5th masses and 2-sigma bounds
%assumes t is (1 x nk), xhatu is (nx x nk), Pu is (nx x nx x nk), x_true is (nx x nk) 
%measurement of position of 5th mass
figs(5)=figure;
ti(5)=tiledlayout(1,2,'TileSpacing','compact','Padding','tight');
nexttile;
plot_estimator(t,xhatu(i1,:),Pu(i1,i1,:),x_true(i1,:),'error');
ylabel('error position q1');
ylim([-0.1 0.1]);
nexttile;
plot_estimator(t,xhatu(i2,:),Pu(i2,i2,:),x_true(i2,:),'error',z);
ylabel('error position q5');
ylim([-0.2 0.2]);
title(ti(5),'(d): constant KF gain with msmt z_{q5}','fontweight','bold');


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