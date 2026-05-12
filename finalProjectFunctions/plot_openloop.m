function plot_openloop(t,x_true,z,x_no_w)
%% plot_openloop(t,x_true,z,x_no_w)
% This function plots a scalar true state along with 
% the measurement (optional) and state with no noise (w=0)
%
%   NOTE: assumes scalar state
% 
% INPUTS
%  t = time vector (1 x n)
%  x_true = vector of the true state (1 x n)
%  z = measurement vector (1 x n) (optional)
%  x_no_w = vector of the state without noise (1 x n) (optional)
% 
%   MAE 6760 Model Based Estimation
%   Cornell University
%   M Campbell
%
if nargin < 4
    x_no_w=[]; %no without noise state given
    if nargin < 3
        z=[]; %no measurement given
    end  
end
%make sure inputs are vectors of the right size
[r1,~]=size(t);[r2,~]=size(x_true);
if (r1~=1) | (r2~=1)
    disp('input vectors are not row vectors')
    return;
end
if abs(t(2)-t(1) - 1) < 1E-10 %t is a vector of integer timesteps
    time_label = 'timestep {\it{k}}';
else
    time_label = 'time {\it{t}} (sec)';
end
%colors for plotting
MCcolors.blue=[4,51,255]/255;
MCcolors.green=[0,160,0]/255;
MCcolors.mag=[255,64,255]/255;
%
pp(1)=plot(t,x_true,'-.','color',MCcolors.mag);
legend_labels={'true state'};
hold on;
%plot measurement
if ~isempty(z)
    pp(2)=plot(t,z,'.','Color',MCcolors.green,'MarkerSize',12);
    legend_labels{end+1}='measurement';
end
%plot no noise case (if used)
if ~isempty(x_no_w)
    pp(3)=plot(t,x_no_w,'-','Color',MCcolors.blue);
    legend_labels{end+1}='state w/ no noise'; 
end
grid;
hold off
%labels and legend
xlabel(time_label);
ylabel('state');
legend(legend_labels,'location','best');

end
