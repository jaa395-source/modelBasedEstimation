function [azDeg, elDeg, range, nedLOS] = ecefRadiusLatLonToAzElRangeNED(obs, tgt)

% Computes azimuth, elevation, range, and NED line-of-sight vector
% from observer to target.
%
% Inputs:
%   rObs       = observer Earth-centered radius
%   latObsDeg  = observer geocentric latitude, deg
%   lonObsDeg  = observer longitude, deg
%   rTgt       = target Earth-centered radius
%   latTgtDeg  = target geocentric latitude, deg
%   lonTgtDeg  = target longitude, deg
%
% Outputs:
%   azDeg   = azimuth, deg, clockwise from north
%   elDeg   = elevation, deg, positive above local horizon
%   range   = straight-line range, same units as radius
%   nedLOS  = line-of-sight vector in local NED frame [N; E; D]

    % Convert both objects to ECEF4
    
    rObs = obs(1);
    lonObsDeg = rad2deg(obs(2));
    latObsDeg = rad2deg(obs(3));

    rTgt = tgt(1);
    lonTgtDeg = rad2deg(tgt(2));
    latTgtDeg = rad2deg(tgt(3));


    rObsECEF = sphToECEF(rObs, latObsDeg, lonObsDeg);
    rTgtECEF = sphToECEF(rTgt, latTgtDeg, lonTgtDeg);

    % Line-of-sight vector from observer to target in ECEF
    losECEF = rTgtECEF - rObsECEF;

    % Local NED basis vectors at observer, expressed in ECEF
    [nHat, eHat, dHat] = nedBasis(latObsDeg, lonObsDeg);

    % Project ECEF LOS into local NED frame
    N = dot(losECEF, nHat);
    E = dot(losECEF, eHat);
    D = dot(losECEF, dHat);

    nedLOS = [N; E; D];

    % Range
    range = norm(nedLOS);

    % Azimuth, clockwise from north
    azDeg = atan2d(E, N);
    azDeg = mod(azDeg, 360.0);
    

    % Elevation, positive upward
    horizontalRange = sqrt(N^2 + E^2);
    elDeg = atan2d(-D, horizontalRange);
end


function rECEF = sphToECEF(radius, latDeg, lonDeg)
% Converts Earth-centered spherical coordinates to ECEF Cartesian position.
%
% Inputs:
%   radius = Earth-centered radius
%   latDeg = geocentric latitude, deg
%   lonDeg = longitude, deg
%
% Output:
%   rECEF = [x; y; z]

    lat = deg2rad(latDeg);
    lon = deg2rad(lonDeg);

    x = radius * cos(lat) * cos(lon);
    y = radius * cos(lat) * sin(lon);
    z = radius * sin(lat);

    rECEF = [x; y; z];
end


function [nHat, eHat, dHat] = nedBasis(latDeg, lonDeg)
% Builds local NED unit vectors at the observer point,
% expressed in ECEF coordinates.

    lat = deg2rad(latDeg);
    lon = deg2rad(lonDeg);

    nHat = [
        -sin(lat) * cos(lon);
        -sin(lat) * sin(lon);
         cos(lat)
    ];

    eHat = [
        -sin(lon);
         cos(lon);
         0.0
    ];

    dHat = [
        -cos(lat) * cos(lon);
        -cos(lat) * sin(lon);
        -sin(lat)
    ];
end