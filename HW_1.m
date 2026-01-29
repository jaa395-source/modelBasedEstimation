%% Problem 1
R_value = 5;
C_value = 0.25;
L_value = 1;

% a
A_matrix = [-R_value/L_value -1/L_value; 1/C_value 0];
B_matrix = [1/L_value 0]';
C_matrix = [0 1];
D_matrix = 0;

% b
T = [0 1; 1/C_value 0];

sys = ss(A_matrix, B_matrix, C_matrix, D_matrix);
[eigenvectors, eigenvalues] = eig(A_matrix)

% c
if rank(obsv(sys)) == size(obsv(sys), 1)
    disp("This system is observable!")
else
    disp("This system is not observable")
end

% d
syms x y
test = [x y];
eig_col_to_remove = 2;
eq1 = test*eigenvectors(:,eig_col_to_remove);
x_value = 1;
y_value = x_value*solve(eq1 == 0)/y;
ic = double([x_value y_value]);

for column = 1:1
    %ic = eigenvectors(:, column);
    [y, tOut, x] = initial(sys, ic');
    figure;
    hold on;
    yyaxis left
    plot(tOut, x(:,1))
    yyaxis right
    plot(tOut, x(:,2));
    for i = 1:6
        xline(pi/abs(min(min(eigenvalues)))*i)
    end
    title("i and v_c vs time with sensor measurment points", 'IC = [' + string(ic(1)) + ", " + string(ic(2)) + "]");
    xlabel("Time (sec)");
    yyaxis left
    ylabel("Current (A)");
    yyaxis right
    ylabel("Voltage (V)");
    hold off;
end
%% Problem 2

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

%2a
C = [1 0 0; 0 1 0]
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
covariance_mat_for_95_percent = cov_mat*4
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
  

%% Problem 3
p_raining = [0.3, 0.7];
p_busy = [0.2, 0.5, 0.3];

yes = p_raining(1);
no = p_raining(2);

l = p_busy(1);
m = p_busy(2);
h = p_busy(3);

combo_running_yes = [0.4 0.2 0.1 0.9 0.7 0.5];
combo_running_no = 1 - combo_running_yes;

% a
probablity_running_any_day = yes*l*combo_running_yes(1) + yes*m*combo_running_yes(2) + ...
    yes*h*combo_running_yes(3) + no*l*combo_running_yes(4) + ...
    no*m*combo_running_yes(5) + no*h*combo_running_yes(6)

% b
probability_running_sunny_day = l*combo_running_yes(4) + ...
    m*combo_running_yes(5) + h*combo_running_yes(6)

probability_running_not_busy_day = yes*combo_running_yes(1) + no*combo_running_yes(4)

if probability_running_sunny_day > probability_running_not_busy_day
    disp("More likely to run when its sunny than when he's not busy, weather is the deciding factor")
else
    disp("More likely to run when he's not bust than when its sunny, schedule is the deciding factor")
end

% c
%P(B,J) = P(J,B)*P(B)/P(A)
probability_of_busy_if_jogging = (yes*combo_running_yes(3) + no*combo_running_yes(6))*h/probablity_running_any_day
