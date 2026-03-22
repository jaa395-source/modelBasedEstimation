%% Example calculating and plotting a 2D covariance ellipse
% 
% MAE 6760 Model Based Estimation
% Cornell University
% M Campbell

P=[6.5 2.5;2.5 6.5];m=[0;0];
[Xe,Ye,U,S,th] = calculate_ellipse(m, P, 1, 50);
plot(Xe,Ye,'b-','linewidth',2);
hold on;
plot(m(1),m(2),'r.','markersize',15);
hold off
legend('1-sigma ellipse','mean');
grid;
title('2D Covariance Ellipse');