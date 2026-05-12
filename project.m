clear; close all; clc;
MAE6760startup; %can adjust font size, figure size at the bottom of this script
global MCcolors;
%% Define Top Level Parameters
load('measurement_data.mat');
tvec = (1:1:size(hypersonicPositionData, 1));
t = tvec;
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

% https://help.agi.com/stk/index.htm#training/DME_Hypersonics.htm?Highlight=Hypersonic

%% Process Az, El, Range Data
% measurementSets = size(allData, 3);
% reformedMeasurements = nan(6, size(allData, 1),  measurementSets);
% 
% for setIdx = 1:measurementSets
%     for timeStep = 1:size(allData, 1)
%         currentLocation = satPositionData(timeStep, :, setIdx);
%         currentLocation = currentLocation(1:3);
%         if timeStep == 1
%             tempData = allData(timeStep, 2:end, setIdx);
%             tempFixedPosition_m = calculate_vehicle_fixed_position(tempData(1), tempData(2), tempData(3), currentLocation);
%         else
%             currData = allData(timeStep, 2:end, setIdx);
%             currFixedPosition_m = calculate_vehicle_fixed_position(currData(1), currData(2), currData(3), currentLocation);
%             currFixedVelocity_m = (currFixedPosition_m - tempFixedPosition_m)/dt;
% 
%             cartesianState = [currFixedPosition_m; currFixedVelocity_m];
% 
%             reformedMeasurements(:, timeStep - 1, setIdx) = return_relevant_state_vector(cartesianState);
%             tempData = currData;
%         end
%     end
% end
% disp("Redefined Data")


for timeStep = 1:size(hypersonicPositionData, 1)
    cartesianState = hypersonicPositionData(timeStep, :);
    x_true(:, timeStep) = return_relevant_state_vector(cartesianState);
end
disp("Redefined Data")

%% Define System
currSet = 1;
nt = size(hypersonicPositionData, 1);
nx = 6;
nz = 4;
nu = 1;
n = nx;

q_scale = 100;
Q = q_scale*[5.0 0.01 0.01 0.01;
             0.01 5.0 0.01 0.01;
             0.01 0.01 5.0 0.01;
             0.01 0.01 0.01 5.0];

Q = q_scale*diag([100, 1, 1, 100]);

R_sigmas = [1000, 0.1, 0.1, 100];
R = diag(R_sigmas);

H = [1 0 0 0 0 0; 
     0 1 0 0 0 0;
     0 0 1 0 0 0;
     0 0 0 1 0 0]; % measure Range, Lat, Lon, and Velocity


v = sqrtm(R)*randn(nz,nt);
z = H*x_true+v;

X0 = x_true(:, 1, currSet);
P0_vals = [10000, 5, 5, 1000, 1, 1];
P0 = diag(P0_vals);

%IF Initial State
Y0 = inv(P0);
y0 = Y0*X0;
yhatu = zeros(nx,nt); yhatu(:,1) = y0;
Yu = zeros(nx,nx,nt); Yu(:,:,1) = Y0;

% define physical state for plotting
xinfo = zeros(nx, nt); xinfo(:,1) = X0;
Pinfo=zeros(nx,nx,nt); Pinfo(:,:,1) = P0;

% invert system matricies once for loop usage
Qi = inv(Q); Ri = inv(R); sigma0 = 0;

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

[f_sym, A_sym, B_sym, X, U, P] = hgv_dynamics_symbolic_no_coriolis();


for k = 1:(nt-1)
    % predict
    currX = xinfo(:,k);
    sigmak = 0;

    A0 = double(subs(A_sym, [X; U; P], [currX; sigmak; params_vec]));
    B0 = double(subs(B_sym, [X; U; P], [currX; sigmak; params_vec]));

    F = diag((ones(1,6))) + A0*dt;
    G = [dt 0 0 0;
         0 dt 0 0;
         0 0 dt 0;
         0 0 0 dt;
         0 0 0 0;
         0 0 0 0];

    Fi = inv(F);
    Mi = Fi'*Yu(:,:,k)*Fi;
    Om = Mi*G*inv(Qi + G'*Mi*G);
    Yp(:,:,k+1) = (eye(nx) - Om*G')*Mi;
    yhatp(:, k+1) = ((eye(nx) - Om*G')*Fi')*yhatu(:,k);



    % update
    yhatu(:,k+1) = yhatp(:,k+1) + H'*Ri*z(:,k+1);
    Yu(:,:,k+1) = Yp(:,:,k+1) + H'*Ri*H;


    % pull out state and covariance for plotting
    Pinfo(:,:,k+1) = inv(Yu(:,:,k+1));
    xinfo(:,k+1) = Pinfo(:,:,k+1)*yhatu(:,k+1);
end
disp("Filter Run")
% Plot
plot_estimator2(t,xinfo,Pinfo,x_true,z,[1 4], {"range (m)", "velocity (m/s)"});
sgtitle('Information Filter (1 msmt)','FontWeight','bold','Fontsize',16);

plot_estimator2(t,xinfo,Pinfo,x_true,z,[2 3], {"latitude (rad)", "longitude (rad)"});
sgtitle('Information Filter (1 msmt)','FontWeight','bold','Fontsize',16);

% plot_estimator2(t,xinfo,Pinfo,x_true,z,[5 6], {"gamma (red)", "psi (rad)"});
% sgtitle('Information Filter (1 msmt)','FontWeight','bold','Fontsize',16);
disp("Data plotted");

%% Plotting Filter Integrity
if plot_integrity
figs(figureCount)=figure('Position',[100 100 800 600]);
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
end


disp("Done!(?)");


%% Define analysis functions
% state vector X = [r, theta, phi, v, gamma, psi]
% ALL ANGLES IN DEGREES


