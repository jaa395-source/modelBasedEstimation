%% 3
p_matrix = [0.8 0.1 0.1;
    0.1 0.6 0.3;
    0.1 0.3 0.6];
meas_number = 2;

% b
ML_function = table(p_matrix(:,meas_number), 'VariableNames',{'z2'}, 'RowNames',{'x1', 'x2', 'x3'})
ML_estimate = string(ML_function.Properties.RowNames(find(ML_function{:,1} == max(ML_function{:,1}))))

%% c
p_x0 = [0.3, 0.3, 0.3];
ML_function_prior = table(ML_function{:,1}.*p_x0', 'VariableNames',{'z2'}, 'RowNames',{'x1', 'x2', 'x3'})
MAP_estimate = string(ML_function_prior.Properties.RowNames(find(ML_function_prior{:,1} == max(ML_function_prior{:,1}))))

%% d
p_x0_2 = [0.85, 0.1, 0.05];
ML_function_new_prior = table(ML_function{:,1}.*p_x0_2', 'VariableNames',{'z2'}, 'RowNames',{'x1', 'x2', 'x3'})
MAP_estimate_2 = string(ML_function_new_prior.Properties.RowNames(find(ML_function_new_prior{:,1} == max(ML_function_new_prior{:,1}))))

%% e
ML_function_two_measurements = table(p_matrix(:,meas_number).*p_matrix(:,meas_number), 'VariableNames',{'z2 x 2'}, 'RowNames',{'x1', 'x2', 'x3'})
ML_estimate = string(ML_function_two_measurements.Properties.RowNames(find(ML_function_two_measurements{:,1} == max(ML_function_two_measurements{:,1}))))

%% f
number_of_measurements = 1;

ML_function_variable_measurements = p_matrix(:,meas_number).^number_of_measurements;
percent_certainties = ML_function_variable_measurements/norm(ML_function_variable_measurements);

while percent_certainties(2) < 0.99

    number_of_measurements = number_of_measurements + 1;
    ML_function_variable_measurements = p_matrix(:,meas_number).^number_of_measurements;
    percent_certainties = ML_function_variable_measurements/norm(ML_function_variable_measurements);

end

answer_string = " It would take " + number_of_measurements + " measurement(s) to be 99% sure the object was a pedestrian";
disp(answer_string);