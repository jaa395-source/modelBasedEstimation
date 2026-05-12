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