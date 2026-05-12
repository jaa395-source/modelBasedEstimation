% Functions for converstion

% Functions to turn measurements into necessary values
function vehicleFixedPosition_m = calculate_vehicle_fixed_position(azimuth, elevation, range, currentLocation)
phi = azimuth;
th = elevation;
r = range;

currentLocation = currentLocation(:);
relCartesian = [r*cosd(th)*cosd(phi);...
    r*cosd(th)*sind(phi);...
    r*sind(th)];

vehicleFixedPosition_m = (currentLocation + relCartesian);
end

function [theta_lat, phi_lon, r_m] = calculate_llr(fixedPosition_m)
earth_radius_m = 6371.8*1000;
lla = ecef2lla(fixedPosition_m(:)');

theta_lat = lla(1);   % degrees
phi_lon = lla(2);   % degrees
r_m = (lla(3) + earth_radius_m);   % m
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

v_ned = rotation_matrix * v(:);

v_north = v_ned(1);
v_east  = v_ned(2);

heading = atan2(v_east, v_north); % radians
psi = mod(rad2deg(heading), 360);

psi = deg2rad(psi);
gamma = deg2rad(gamma);

v_speed = v_norm;

end

function hypersonicStateVector = return_relevant_state_vector(cartesianStateVector)

currFixedPosition_m = cartesianStateVector(1:3);
[theta_lat, phi_lon, r_m] = calculate_llr(currFixedPosition_m);
[v_speed, gamma, psi] = calculate_velocity_values(cartesianStateVector);

hypersonicStateVector = [r_m; deg2rad(theta_lat); deg2rad(phi_lon); v_speed; gamma; psi];

end