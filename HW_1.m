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
[eigenvalues, eigenvectors] = eig(A_matrix)

% c
if rank(obsv(sys)) == size(obsv(sys), 1)
    disp("This system is observable!")
else
    disp("This system is not observable")
end

% d
for column = 1:2
    ic = eigenvalues(:, column);
    [y, tOut, x] = initial(sys, ic);
    figure;
    hold on;
    yyaxis left
    plot(tOut, x(:,1))
    yyaxis right
    plot(tOut, x(:,2));
    for i = 1:6
        xline(pi/abs(min(min(eigenvectors)))*i)
    end
    title("i and v_c vs time with sensor measurment points", 'IC = [' + string(ic(1)) + ", " + string(ic(2)) + "]");
    xlabel("Time (sec)");
    yyaxis left
    ylabel("Current (A)");
    yyaxis right
    ylabel("Voltage (V)");
    hold off;
end


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
