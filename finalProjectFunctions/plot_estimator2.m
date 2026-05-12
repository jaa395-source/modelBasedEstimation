function fig = plot_estimator2(t,xhat,Phat,x_true,z,ii, names, axis_lims);
%
fig = figure('Position',[100 100 1600 600]);
tiledlayout(1,2,'TileSpacing','compact','Padding','tight');
nexttile;
plot_estimator(t,xhat(ii(1),:),Phat(ii(1),ii(1),:),x_true(ii(1),:),'error',z(ii(1),:))
ylabel(string(names(1)));
axis(axis_lims);
nexttile;
plot_estimator(t,xhat(ii(2),:),Phat(ii(2),ii(2),:),x_true(ii(2),:),'error',z(ii(2),:))
ylabel(string(names(2)));
axis(axis_lims);
%
end