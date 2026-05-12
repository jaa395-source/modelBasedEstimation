clc; clear; close all;

t = linspace(0, 2*pi, 100);
scale = sqrt(2);
sigma_x = 10;
sigma_y = 5;

x_c = 0;
y_c = 0;

x = x_c + (sigma_x)*cos(t);
y = y_c + (sigma_y)*sin(t);

figure;
hold on;
plot(x,y);


x = x_c + (sigma_x*2)*cos(t);
y = y_c + (sigma_y*2)*sin(t);
plot(x,y);
axis equal