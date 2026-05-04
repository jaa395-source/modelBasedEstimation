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
