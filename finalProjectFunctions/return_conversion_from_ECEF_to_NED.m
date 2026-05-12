function R = return_conversion_from_ECEF_to_NED(lat, lon)
R = [-sin(lat)*cos(lon), -sin(lat)*sin(lon),  cos(lat);
    -sin(lon),           cos(lon),          0;
    -cos(lat)*cos(lon), -cos(lat)*sin(lon), -sin(lat)];
end