function [theta_lat, phi_lon, r_m] = calculate_llr(fixedPosition_m)
earth_radius_m = 6371.8*1000;
lla = ecef2lla(fixedPosition_m(:)');

theta_lat = lla(1);   % degrees
phi_lon = lla(2);   % degrees
r_m = (lla(3) + earth_radius_m);   % m
end