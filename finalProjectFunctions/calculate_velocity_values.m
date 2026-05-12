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