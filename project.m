%% Define Top Level Parameters
hypersonicVehicleTrajectoryFilePath_J2000 = "filename";
drag_coeff = 1.0;
lift_coeff = 1.0;
surface_area = 1.0;
mass = 1.0;
C_gamma = 0;
C_psi = 0;
noise  = 0;
dt = 1;



%% Load in data
trueTrajectory = readtable(hypersonicVehicleTrajectoryFilePath_J2000);
initPosition = trueTrajectory(1,:);



% Plot true trajectory
figure();
hold on;
title("True vs Calculated Trajectory");
predicted_steps = size(trueTrajectory, 1) - 1;

nextSteps = initPosition;
for predIdx = 1:predicted_steps

    currentStep = nextSteps(end, :);
    currentLocation = currentStep(1:3);
    currentVelocity = currentStep(4:6);

    rangeToEarth_km = norm(currentLocation);
    dynamicPressure = calculate_dynamic_pressure(rangeToEarth_km, currentVelocity);
    accelerations = calculate_flight_accelerations(dynamicPressure, drag_coeff, lift_coeff, surface_area, mass);

    propagatedStep = propagate_state(currentStep, noise,...
        accelerations, sigma_deg, C_gamma, C_psi, dt);
    nextSteps(size(nextSteps, 1) + 1, :) = propagatedStep;
    




end
%% Define analysis functions
% state vector X = [r, theta, phi, v, gamma, psi]
% ALL ANGLES IN DEGREES
function rangeToEarth_km = calculate_range_to_earth(azimuth, elevation, range, currentLocation)
phi = azimuth;
th = elevation;
r = range;

currentLocation = currentLocation(:);
relCartesian = [r*cosd(th)*cosd(phi);...
                r*cosd(th)*cosd(phi);...
                r*sind(th0)];

rangeToEarth_km = norm(currentLocation + relCartesian);
end


function dynamicPressure = calculate_dynamic_pressure(rangeToEarth_km, v)
rho_0 = 1.225; %kg/m^3
earth_radius_km = 6371.8;
h_m = (rangeToEarth_km - earth_radius_km)*1000;
hs_m = 6700;
rho = rho_0.*exp(-h_m/hs_m);

dynamicPressure = rho.*(v.*v)/2;
end

function accelerations = calculate_flight_accelerations(dynamicPressure, drag_coeff, lift_coeff, surface_area, mass)
    accelerations = [dynamicPressure.*drag_coeff.*surface_area./mass;
                     dynamicPressure.*lift_coeff.*surface_area./mass];
end

function Xk = propagate_state(Xk, noise, accelerations, sigma_deg, C_gamma, C_psi, dt)
r = Xk(1,:);
th = Xk(2,:);
phi = Xk(3,:);
v = Xk(4,:);
gamma = Xk(5,:);
psi = Xk(6,:);

drag = accelerations(1,:);
lift = accelerations(2,:);


r_dot = v.*sind(gamma);
th_dot = v.cosd(gamma).*sind(psi)./r./cosd(phi);
phi_dot = v.cosd(gamma).*cosd(psi)./r;
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