clc; clear; close all;
%
%   Homework #1
%   problem #2: Multivariate Gaussian Random Variable
%        application: localizat; clear; close all; 
%% MAE ion of a phone
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
convert_LLerr2merr = (pi/180)*Re; 

C = [1 0 0; 0 1 0];
x_ll = C*x_lla;
x_ll_temp = [0 0]';

cov_mat = C*P_llam;
cov_mat = cov_mat(:,1:2);

% 2b
[Xe,Ye,U,S,th] = calculateEllipseCovMC(x_ll_temp, cov_mat, 1, 36);
figure;
hold on;
plot(Xe, Ye);
ylabel("East/West Error (m)");
xlabel("North/South Error (m)");
plot(0,0, '+', 'LineWidth',4)

max_xe = max(Xe);
pair_ye = Ye(find(Xe == max(Xe)));
semi_major_axis_68 = (norm([max_xe pair_ye]))

max_ye = max(Xe);
pair_xe = Xe(find(Ye == max(Ye)));
semi_minor_axis_68 = norm([max_ye pair_xe])
legend("68% Confidence Elipse", "Center Point");

%2c
covariance_mat_for_95_percent = cov_mat*4;
[Xe,Ye,U,S,th] = calculateEllipseCovMC(x_ll_temp, covariance_mat_for_95_percent, 1, 36);
plot(Xe, Ye);
max_xe = max(Xe);
pair_ye = Ye(find(Xe == max(Xe)));
semi_major_axis_95 = norm([max_xe pair_ye])

max_ye = max(Xe);
pair_xe = Xe(find(Ye == max(Ye)));
semi_minor_axis_95 = norm([max_ye pair_xe])

axis equal
legend("68% Confidence Elipse", "Center Point", "95% Confidence Elipse");
% maintain the same defintion for error and use the new covariance matrix,
% which 4 times the original 1 sigma covariance matrix

%2d
e_m = (x_true - x_ll)*convert_LLerr2merr;
f_value = e_m'*inv(covariance_mat_for_95_percent)*e_m;
if f_value < 1
    answer_string = "Estimator is working well, f < 1";
else
    answer_string = "Estimator is NOT working well, f > 1";
end
disp(answer_string);
  