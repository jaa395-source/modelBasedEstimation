%% Descriptions
% measurements - azimuth, elevation, range, speed
% z = 4x1
% x = 6x1
clear all; close all; clc;
rng(10);
%% Load in Functions
addpath("finalProjectFunctions\");

%% Load in true state and define parameters
load('measurement_data.mat');
tvec = (1:1:size(hypersonicPositionData, 1));
nt = length(tvec);
plot_integrity = false;

params.Re = 6371.1*1000;
params.hs = 6700;
params.rho0 = 1.225;
params.CD= 0.015;
params.CL = 0.085;
params.S = 3.99;
params.m = 1400.0;
params.mu = 3.99e14;
params.Cgamma = 0;
params.Cpsi = 0;
params.omegaE = 7.2921159e-5;   % Earth rotation rate, rad/s
noise  = 0;
dt = 0.99;

params_vec = [
    6371.1e3;
    3.986004418e14;
    1.225;
    6700;
    3.99;
    1400;
    0.085;
    0.015
    ];
disp("Defined Parameters");

%% Prepare data for analysis
measurementSet = 1;
for timeStep = 1:nt
    cartesianState = hypersonicPositionData(timeStep, :);
    x_true(:, timeStep) = return_relevant_state_vector(cartesianState);
end
disp("Generated simulated state");

z_clean = x_true;
nz = 4;
H = [1 0 0 0 0 0;
    0 1 0 0 0 0;
    0 0 1 0 0 0;
    0 0 0 1 0 0];
H = [H zeros(4,2)];

q_scale = 100;
Q = q_scale*diag([50, 50 0.1 0.1]);

R_sigmas = [1000, 0.5, 0.5, 100];
R = diag(R_sigmas);

v = sqrtm(R)*randn(nz,nt);
long_bias = 7.29e-5;
lat_bias = 0.005;

x_true=[x_true;long_bias*ones(1,nt);lat_bias*ones(1,nt)];
z = (x_true(1:4,:)) + v;

nw = 4;
w = sqrtm(Q)*randn(nw,nt);
drag_values = w(1,:);
lift_values = w(2,:);
long_bias_values = w(3,:);
lat_bias_values = w(4,:);
[f_sym, A_sym, B_sym, X, U, P] = hgv_dynamics();

%% Define relevant values
nk = size(x_true,2);
nx = 6;
nz = 4;
nw = 2;

P0_vals = [];
for i = 1:nx
    P0_vals(i) = cov(x_true(i,:));
end

P0 = diag([P0_vals, 0.5, 0.5]);
nx=nx+2;

% Kalman
x0 = [x_true(:,1)];
xhatp = zeros(nx,nk);xhatp(:,1) = x0;
xhatu = zeros(nx,nk);xhatu(:,1) = x0;
Pp = zeros(nx,nx,nk);Pp(:,:,1) = P0;
Pu = zeros(nx,nx,nk);Pu(:,:,1) = P0;
n = nk;

% Information

%IF Initial State
Y0 = inv(P0);
y0 = Y0*x0;
yhatu = zeros(nx,nk); yhatu(:,1) = y0;
Yu = zeros(nx,nx,nk); Yu(:,:,1) = Y0;

% define physical state for plotting
xinfo = zeros(nx, nk); xinfo(:,1) = x0;
Pinfo=zeros(nx,nx,nk); Pinfo(:,:,1) = P0;

% invert system matricies once for loop usage
Qi = inv(Q); Ri = inv(R);

% msmt gating parameters
Pgate = 0.95;alpha = 1-Pgate;
Lam0  = chi2inv(Pgate, nz);

%filter consistency parameters
PF = 0.95;alpha = 1-PF;
win=10;%time window
Blow=chi2inv(alpha/2,10*nz)/win;  %low filter threshold
Bhigh=chi2inv(1-alpha/2,10*nz)/win;  %high filter threshold
Lam=zeros(nt,1)';
LamF=zeros(nt,1)';
NFrej=0; TFrej=[]; EFrej=[];
Nrej=0;  Trej=[];  Erej=[];


for k=1:(nk-1)

    % Kalman predict
    currX = xhatu(:,k);
    uk = [drag_values(k); lift_values(k); long_bias_values(k); lat_bias_values(k)];
    A0 = double(subs(A_sym, [X; U; P], [currX; uk; params_vec]));
    B0 = double(subs(B_sym, [X; U; P], [currX; uk; params_vec]));

    [F,G]=c2d(A0,B0,dt);

    xhatp(:,k+1) = F*currX + G*uk;
    Pp(:, :, k+1) = F*Pu(:,:,k)*F' + G*Q*G';

    % Kalman gain
    K = Pp(:,:,k+1)*H'*inv(H*Pp(:,:,k+1)*H' + R);

    % Kalman Update
    xhatu(:,k+1) = xhatp(:,k+1) + K*(z(:,k+1) - H*xhatp(:,k+1));
    Pu(:,:,k+1) = (eye(nx) - K*H)*Pp(:,:,k+1)*(eye(nx) - K*H)' + K*R*K';


    % Information predict
    Fi = inv(F);
    Mi = Fi'*Yu(:,:,k)*Fi;
    Om = Mi*G*inv(Qi + G'*Mi*G);
    Yp(:,:,k+1) = (eye(nx) - Om*G')*Mi;
    yhatp(:, k+1) = ((eye(nx) - Om*G')*Fi')*yhatu(:,k);

    % Information update
    yhatu(:,k+1) = yhatp(:,k+1) + H'*Ri*z(:,k+1);
    Yu(:,:,k+1) = Yp(:,:,k+1) + H'*Ri*H;

    % Information pull out state and covariance for plotting
    Pinfo(:,:,k+1) = inv(Yu(:,:,k+1));
    xinfo(:,k+1) = Pinfo(:,:,k+1)*yhatu(:,k+1);

    % Consistency Calculations
    % innovation
    inn = (z(:, k+1) - H*xhatp(:, k+1));
    S = H*Pp(:, :, k + 1)*H' + R;

    Lam(1,k+1) = inn'*inv(S)*inn;

    if Lam(1,k+1)>Lam0
        Nrej = Nrej + 1;
        Trej(length(Trej) + 1) = tvec(k+1);
        Erej(:,size(Erej,2) + 1) = z(:,k+1)-H*x_true(:,k+1);
    end

    if k>=win
        LamF(1,k+1) = mean(Lam(k-win+2:k+1));
    end
    if (LamF(1,k+1) < Blow) | (LamF(k+1)>Bhigh)
        NFrej = NFrej + 1;
        TFrej(length(TFrej) + 1) = tvec(k+1);
        EFrej(:,size(EFrej,2) + 1) = z(:,k+1)-H*x_true(:,k+1);
    end
end

%% Plot Initial Comparison
axis_lims = [0 440 -75 75];
plot_estimator2(tvec,xhatu,Pu,x_true,z,[1 4], {"range (m)", "velocity (m/s)"}, axis_lims);
sgtitle('Kalman Filter','FontWeight','bold','Fontsize',16);

axis_lims = [0 440 -2 2];
plot_estimator2(tvec,xhatu,Pu,x_true,z,[2 3], {"longitude (rad)", "latitude (rad)"}, axis_lims);
sgtitle('Kalman Filter','FontWeight','bold','Fontsize',16);

axis_lims = [0 440 -75 74];
plot_estimator2(tvec,xinfo,Pinfo,x_true,z,[1 4], {"range (m)", "velocity (m/s)"}, axis_lims);
sgtitle('Information Filter (1 msmt)','FontWeight','bold','Fontsize',16);

axis_lims = [0 440 -2 2];
plot_estimator2(tvec,xinfo,Pinfo,x_true,z,[2 3], {"longitude (rad)", "latitude (rad)"}, axis_lims);
sgtitle('Information Filter (1 msmt)','FontWeight','bold','Fontsize',16);


%% Run Information Filter with multiple sensors
R_sigmas = [1000, 0.5, 0.5, 100];
R = diag(R_sigmas);

for setIdx = 2:size(allData, 3)
    v = sqrtm(R)*randn(nz,nt);
    allZ(:,:,setIdx - 1) = (x_true(1:4,:)) + v;
end

q_scale = 100;
Q = q_scale*diag([50, 50 0.1 0.1]);
Q = 1.7*Q;
%IF Initial State
Y0 = inv(P0);
y0 = Y0*x0;
yhatu = zeros(nx,nk); yhatu(:,1) = y0;
Yu = zeros(nx,nx,nk); Yu(:,:,1) = Y0;

% define physical state for plotting
xinfo_many = zeros(nx, nk); xinfo_many(:,1) = x0;
Pinfo_many=zeros(nx,nx,nk); Pinfo_many(:,:,1) = P0;

% invert system matricies once for loop usage
Qi = inv(Q); Ri = inv(R);

% msmt gating parameters
Pgate = 0.95;alpha = 1-Pgate;
Lam0  = chi2inv(Pgate, nz);

%filter consistency parameters
PF = 0.95;alpha = 1-PF;
win=10;%time window
Blow=chi2inv(alpha/2,10*nz)/win;  %low filter threshold
Bhigh=chi2inv(1-alpha/2,10*nz)/win;  %high filter threshold
Lam=zeros(nt,1)';
LamF=zeros(nt,1)';
NFrej=0; TFrej=[]; EFrej=[];
Nrej=0;  Trej=[];  Erej=[];

for k=1:(nk-1)

    % Kalman predict
    currX = xhatu(:,k);
    uk = [drag_values(k); lift_values(k); long_bias_values(k); lat_bias_values(k)];
    A0 = double(subs(A_sym, [X; U; P], [currX; uk; params_vec]));
    B0 = double(subs(B_sym, [X; U; P], [currX; uk; params_vec]));

    [F,G]=c2d(A0,B0,dt);

    % Information predict
    Fi = inv(F);
    Mi = Fi'*Yu(:,:,k)*Fi;
    Om = Mi*G*inv(Qi + G'*Mi*G);
    Yp(:,:,k+1) = (eye(nx) - Om*G')*Mi;
    yhatp(:, k+1) = ((eye(nx) - Om*G')*Fi')*yhatu(:,k);

    % Information update
    yhatu(:,k+1) = yhatp(:,k+1) + H'*Ri*z(:,k+1);
    Yu(:,:,k+1) = Yp(:,:,k+1) + H'*Ri*H;

    for setIdx = 1:size(allZ, 3)
        yhatu(:,k+1) = yhatu(:,k+1) + H'*Ri*allZ(:,k+1, setIdx);
        Yu(:,:,k+1) = Yu(:,:,k+1) + H'*Ri*H;
    end

    % Information pull out state and covariance for plotting
    Pinfo_many(:,:,k+1) = inv(Yu(:,:,k+1));
    xinfo_many(:,k+1) = Pinfo_many(:,:,k+1)*yhatu(:,k+1);

    % Consistency Calculations
    % innovation
    inn = (z(:, k+1) - H*xhatp(:, k+1));
    P1p = inv(Yp(:,:,k+1));
    S = H*P1p*H' + R;

    Lam(1,k+1) = inn'*inv(S)*inn;

    if Lam(1,k+1)>Lam0
        Nrej = Nrej + 1;
        Trej(length(Trej) + 1) = tvec(k+1);
        Erej(:,size(Erej,2) + 1) = z(:,k+1)-H*x_true(:,k+1);
    end

    if k>=win
        LamF(1,k+1) = mean(Lam(k-win+2:k+1));
    end
    if (LamF(1,k+1) < Blow) | (LamF(k+1)>Bhigh)
        NFrej = NFrej + 1;
        TFrej(length(TFrej) + 1) = tvec(k+1);
        EFrej(:,size(EFrej,2) + 1) = z(:,k+1)-H*x_true(:,k+1);
    end
end

msmt_count = string(size(allData, 3));
axis_lims = [0 440 -30 30];
plot_estimator2(tvec,xinfo_many,Pinfo_many,x_true,z,[1 4], {"range (m)", "velocity (m/s)"}, axis_lims);
sgtitle('Information Filter (' + msmt_count + ' msmts)','FontWeight','bold','Fontsize',16);

axis_lims = [0 440 -0.7 0.7];
plot_estimator2(tvec,xinfo_many,Pinfo_many,x_true,z,[2 3], {"longitude (rad)", "latitude (rad)"}, axis_lims);
sgtitle('Information Filter (' + msmt_count + ' msmts)','FontWeight','bold','Fontsize',16);

disp("Information filter run for " + msmt_count + " sensors");

%
%uncomment lines below to output the percent of msmts rejected
disp('(c) percent of msmts rejected:');
disp(Nrej/nk*100)
%uncomment lines below to output the percent of time filter is inconsistent
disp('(c) percent of time filter is inconsistent:');
disp(NFrej/nk*100)



