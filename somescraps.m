% Scraps
%% Load in data
trueTrajectory = readtable(hypersonicVehicleTrajectoryFilePath_J2000);
firstUsedDataPoint = 1;
totalDataPoints = size(trueTrajectory, 1);
totalUsedDataPoints = length(firstUsedDataPoint:totalDataPoints);
initPosition = trueTrajectory{firstUsedDataPoint,:}*1000;

        [theta_lat, phi_lon, r_alt_m] = calculate_llr(initPosition(2:4));
        [v_speed, gamma, psi] = calculate_velocity_values(initPosition(2:end));
        X0 = [r_alt_m; theta_lat; phi_lon; v_speed; gamma; psi];



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