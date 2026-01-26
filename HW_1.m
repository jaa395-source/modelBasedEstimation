%% Problem 1
R_value = 5;
C_value = 0.25;
L_value = 1;

% a
A_matrix = [-1/L_value -R_value/L_value; 1/L_value 0];
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
column = 1;
ic = eigenvalues(:, column);
[y, tOut, x] = initial(sys, ic);
