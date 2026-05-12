function hypersonicStateVector = return_relevant_state_vector(cartesianStateVector)

currFixedPosition_m = cartesianStateVector(1:3);
[theta_lat, phi_lon, r_m] = calculate_llr(currFixedPosition_m);
[v_speed, gamma, psi] = calculate_velocity_values(cartesianStateVector);

hypersonicStateVector = [r_m; deg2rad(theta_lat); deg2rad(phi_lon); v_speed; gamma; psi];

end