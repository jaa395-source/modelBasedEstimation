clear; close all; clc;
%% Define Top Level Parameters
hypersonicVehicleTrajectoryFilePath_J2000 = "HXRV_X43_PosVel.csv";
drag_coeff = 0.015;
lift_coeff = 0.085;
surface_area = 3.99;
mass = 1400.0;
C_gamma = 0;
C_psi = 0;
noise  = 0;
dt = 1;

% https://help.agi.com/stk/index.htm#training/DME_Hypersonics.htm?Highlight=Hypersonic

%% Process Az, El, Range Data
measurementSets = size(allData, 3);
reformedMeasurements = nan(6, size(allData, 1) - 1,  measurementSets);

for setIdx = 1:measurementSets
    for timeStep = 1:size(allData, 1)
        currentLocation = positionData(timeStep, :, setIdx);
        %currentLocation = [0;0;0];
        if timeStep == 1
            tempData = allData(timeStep, 2:end, setIdx);
            tempFixedPosition_m = calculate_vehicle_fixed_position(tempData(1), tempData(2), tempData(3), currentLocation);
        else
            currData = allData(timeStep, 2:end, setIdx);
            currFixedPosition_m = calculate_vehicle_fixed_position(currData(1), currData(2), currData(3), currentLocation);
            currFixedVelocity_m = (currFixedPosition_m - tempFixedPosition_m)/dt;
            reformedMeasurements(:, timeStep, measurementSets) = [currFixedPosition_m; currFixedVelocity_m];
            tempData = currData;
        end
    end
end

X_true = zeros(6, totalUsedDataPoints);
X_pred = zeros(6, totalUsedDataPoints);

X_true(:,1) = X0;
X_pred(:,1) = X0;

for predIdx = firstUsedDataPoint+1:predicted_steps
    shiftedIdx = predIdx - firstUsedDataPoint + 1;

    % Recalculate true values
    trueIdx = trueTrajectory{predIdx,:}*1000;
    [theta_lat, phi_lon, r_alt_m] = calculate_llr(trueIdx(2:4));
    [v_speed, gamma, psi] = calculate_velocity_values(trueIdx(2:end));
    X_true(:, shiftedIdx) = [r_alt_m; theta_lat; phi_lon; v_speed; gamma; psi];

    % Calculate predicted values
    if use_predicition
        X_prev = X_pred(:, shiftedIdx - 1);
    else
        X_prev = X_true(:, shiftedIdx);
    end
    rangeToEarth_m = X_prev(1);
    currentVelocity_m_sec = X_prev(4);
    dynamicPressure = calculate_dynamic_pressure(rangeToEarth_m, currentVelocity_m_sec);
    accelerations = calculate_flight_accelerations(dynamicPressure, drag_coeff, lift_coeff, surface_area, mass);

    propagatedStep = propagate_state(X_prev, noise,...
        accelerations, C_gamma, C_psi, dt);
    X_pred(:, shiftedIdx) = propagatedStep;

end


% Plot true trajectory
figure();
hold on;
title("True vs Calculated Trajectory");
predicted_steps = size(trueTrajectory, 1) - 1;
use_predicition = true;
plot(X_true(3,:), X_true(2,:));
plot(X_pred(3,:), X_pred(2,:));
legend({'True', 'Predicted'});
hold off;


% Calculate Residuals
figure();
hold on;
title("Residual Ratios vs Time");
xlabel("Time (EpSec)");
ylabel("Residual Ratios");

residuals = X_true - X_pred;
lat_residual_ratios = residuals(2,:)./X_true(2,:);
lon_residual_ratios = residuals(3,:)./X_true(3,:);
plot(lat_residual_ratios);
plot(lon_residual_ratios);
legend({'Latitude Residual Ratios', 'Longitude Residual Ratios'});

disp("Done!(?)");
%% Define analysis functions
% state vector X = [r, theta, phi, v, gamma, psi]
% ALL ANGLES IN DEGREES

% Functions for converstion
function R = return_conversion_from_ECEF_to_NED(lat, lon)
R = [-sin(lat)*cos(lon), -sin(lat)*sin(lon),  cos(lat);
    -sin(lon),           cos(lon),          0;
    -cos(lat)*cos(lon), -cos(lat)*sin(lon), -sin(lat)];
end

% Functions to turn measurements into necessary values
function vehicleFixedPosition_m = calculate_vehicle_fixed_position(azimuth, elevation, range, currentLocation)
phi = azimuth;
th = elevation;
r = range;

currentLocation = currentLocation(:);
relCartesian = [r*cosd(th)*cosd(phi);...
    r*cosd(th)*cosd(phi);...
    r*sind(th)];

vehicleFixedPosition_m = (currentLocation + relCartesian);
end

function [theta_lat, phi_lon, r_m] = calculate_llr(fixedPosition_m)
earth_radius_m = 6371.8*1000;
lla = ecef2lla(fixedPosition_m);

theta_lat = lla(1);   % degrees
phi_lon = lla(2);   % degrees
r_m = (lla(3) + earth_radius_m);   % km
end

function [v_speed, gamma, psi] = calculate_velocity_values(fixedState)

r_m = fixedState(1:3);        % position (m)
v = fixedState(4:6);     % velocity (m/s)

r_norm = norm(r_m);
v_norm = norm(v);

gamma_rad = asin(dot(r_m, v) / (r_norm * v_norm)); % radians
gamma = rad2deg(gamma_rad);

[theta_lat, phi_lon, ~] = calculate_llr(r_m);
rotation_matrix = return_conversion_from_ECEF_to_NED(deg2rad(theta_lat), deg2rad(phi_lon));

v_ned = rotation_matrix * v';

v_north = v_ned(1);
v_east  = v_ned(2);

heading = atan2(v_east, v_north); % radians
psi = mod(rad2deg(heading), 360);

v_speed = v_norm;

end

% Convert Az, El, Range to state vector
function Xk = calculate_state_vector_from_measurements(fixedState)
fixedPosition_m = fixedState(1:3);

[theta_lat, phi_lon, r_m] = calculate_llr(fixedPosition_m);
[v_speed, gamma, psi] = calculate_velocity_values(fixedState);

Xk = [r_m; theta_lat; phi_lon; v_speed; gamma; psi];
end


% Functions to calculate propogation values
function dynamicPressure = calculate_dynamic_pressure(rangeToEarth_m, v)
rho_0 = 1.225; %kg/m^3
earth_radius_m = 6371.8*1000;
h_m = (rangeToEarth_m - earth_radius_m);
hs_m = 6700;
rho = rho_0.*exp(-h_m/hs_m);

dynamicPressure = rho.*(v.*v)/2;
end

function accelerations = calculate_flight_accelerations(dynamicPressure, drag_coeff, lift_coeff, surface_area, mass)
accelerations = [dynamicPressure.*drag_coeff.*surface_area./mass;
    dynamicPressure.*lift_coeff.*surface_area./mass];
end

% Propogation Functions
function Xk = propagate_state(Xk, noise, accelerations, C_gamma, C_psi, dt)
sigma_deg = 35;
g = 9.81;
r = Xk(1,:);
th = Xk(2,:);
phi = Xk(3,:);
v = Xk(4,:);
gamma = Xk(5,:);
psi = Xk(6,:);

drag = accelerations(1,:);
lift = accelerations(2,:);

r_dot = v.*sind(gamma);
th_dot = v.*cosd(gamma).*sind(psi)./r./cosd(phi);
phi_dot = v.*cosd(gamma).*cosd(psi)./r;
v_dot = -drag - g.*sind(gamma);
gamma_dot = (1./v).*(lift.*cosd(sigma_deg) +  (((v.*v)./r) - g).*cosd(gamma)) + C_gamma;
psi_dot = (1./v).*((lift.*sind(sigma_deg)./cosd(gamma)) +  ((v.*v)./r).*cosd(gamma).*sind(psi).*tand(phi)) + C_psi;

Xk = Xk + dt*[r_dot;
    th_dot;
    phi_dot;
    v_dot;
    gamma_dot;
    psi_dot];
end