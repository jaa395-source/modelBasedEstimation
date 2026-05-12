%% MAE 6760 Model Based Estimation
%
%   Homework #1
%   problem #4: Linear System driven by Random Processes
%        application: quarter suspension over a bumpy road
%

%constants
M=300;    %car mass (kg)
m=50;     %wheel mass (kg)
K1=3000;  %spring constant (N/m)
K2=30000; %spring constant (N/m)
C1=600;   %damping constant (Nsec/m)

% Set up open loop model, including r,u as inputs and z as output
%
% xdot = A*x + Bu*u + Br*r
%    z = C*x + Du*u + Dr*r
%

%State Space system matrices
A=[0  0  1  0; 
   0  0  0  1 ; 
   -K1/M K1/M -C1/M C1/M ;
   K1/m -(K1+K2)/m C1/m -C1/m];
Bu=[0; 0; 1/M; -1/m];  %The actuator control as an input u(t)
Br=[0; 0; 0; K2/m];    %The bumpy road as an input r(t)
C=[1 0 0 0];Du=0;Dr=0; %The driver position as the output
% Bumpy Road white noise disturbance intensity
Sigr = 2E-4; %m^2/sec

%% -------------------
%% Part (a): simulate open loop system
Tf=1000;dt=0.01;t=[0:dt:Tf]';
r=sqrt(Sigr)*randn(length(t),1)/sqrt(dt);

sys_CT = ss(A, Br, C, Dr);
[ol_z_CT, ~, ol_x_CT] = lsim(sys_CT, r, t);
ol_Pz_CT = cov(ol_z_CT)

Px_CTan = lyap(A, Br*Sigr*Br');
ol_Pz_CTan = C*Px_CTan*C'

percent_difference = ((ol_Pz_CTan - ol_Pz_CT)/(ol_Pz_CT))*100

%% -------------------
%% Part (b): find and simulate closed loop system
% These lines find the closed loop state feedback controller K,
% where the form of the controller is: u = -K*x
Rzz=1;Ruu=2E-9; 
[K,S,E]=lqry(ss(A,Bu,C,Du),Rzz,Ruu);
%  
%Simulate the closed loop system for the bumpy road
cl_sys = ss(A-Bu*K, Br, C, Dr);
[cl_z_CT, ~, cl_x_CT] = lsim(cl_sys, r, t);
cl_Pz_CT = cov(cl_z_CT)

Px_CTan = lyap(A-Bu*K, Br*Sigr*Br');
cl_Pz_CTan = C*Px_CTan*C'

cl_percent_difference = ((cl_Pz_CTan - cl_Pz_CT)/(cl_Pz_CT))*100


%% -------------------
%% Part (c): plot open and closed loop response z(t), analyze
figure;
hold on;
plot(t, ol_z_CT, t, cl_z_CT);
title("Z value vs Time");
xlabel("Time (sec)");
ylabel("Performance");
legend("Open Loop Data", "Closed Loop Data");

ol_stdv = std(ol_z_CT);
cl_stdv = std(cl_z_CT);

limit_meters = .03;
p = normcdf([-limit_meters/ol_stdv limit_meters/ol_stdv]);
probablity_of_more_than_3_cm_ol_percent = (1 - (p(2) - p(1)))*100

p = normcdf([-limit_meters/cl_stdv limit_meters/cl_stdv]);
probablity_of_more_than_3_cm_cl_percent = (1 - (p(2) - p(1)))*100



%% -------------------
%% Part (d): analyze the control effort u(t)


