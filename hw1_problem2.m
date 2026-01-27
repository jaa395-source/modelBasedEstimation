%% MAE 6760 Model Based Estimation
%
%   Homework #1
%   problem #2: Multivariate Gaussian Random Variable
%        application: localization of a phone
%

%Upson 206 "true" Lat-Lon location
x_true=[ 42.44396; %deg
        -76.48248];%deg

Re=6378100; %radius of the Earth in m

%3D position in LLA coordinates
x_lla=[42.44445; %deg
      -76.48252; %deg
            263];%m
%3D covariance position in units of m^2
P_llam = [400 40 100;
          40 400 100;
          100 100 2500]; 

%conversion for Lat-Lon angular error to error in m
%assumes small angles
convert_LLerr2merr = pi/180*Re; 