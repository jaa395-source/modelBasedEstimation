%% Define Top Level Parameters
hypersonicVehicleTrajectoryFilePath_J2000 = "filename";


%% Define analysis functions
% state vector X = [r, theta, phi, v, gamma, psi]
% ALL ANGLES IN DEGREES
function rangeToEarth = calculate_range_to_earth(azimuth, elevation, range, currentLocation)
phi = azimuth;
th = elevation;
r = range;

currentLocation = currentLocation(:);
relCartesian = [r*cosd(th)*cosd(phi);...
                r*cosd(th)*cosd(phi);...
                r*sind(th0)];

rangeToEarth = norm(currentLocation + relCartesian);
end

function Xk = propagate_state(Xk, noise, drag, lift, sigma_deg, C_gamma, C_psi, dt)

r = Xk(1,:);
th = Xk(2,:);
phi = Xk(3,:);
v = Xk(4,:);
gamma = Xk(5,:);
psi = Xk(6,:);


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