function plot_estimator(t,xest,Pest,x_true,plot_type,z)
%% plot_estimator(t,xest,Pest,x_true,plot_type,z)
% This function plots a scalar state estimate and the +/- 2-sigma bounds 
% from a filter along with the true value (optional) and measurement (optional)
% 
%   NOTE: assumes scalar state
% 
% INPUTS
%  t = time vector (1 x n)
%  xest = vector of state estimates (1 x n)
%  Pest = vector of error variance of estimator (1 x n)
%  x_true = vector of true state (1 x n)
%  plot_type = 'state':     plots xest +/- 2*sigma, x_true
%              'error':     plots err +/- 2*sigma
%  z = measurement vector (1 x n) (optional)
% 
%   MAE 6760 Model Based Estimation
%   Cornell University
%   M Campbell
%
if nargin < 6
    z=[]; %no measurement given
    if nargin < 5
        plot_type = 'error'; %default
    end  
end
%make sure inputs are vectors of the right size
[r1,~]=size(t);[r2,~]=size(xest);[r3,~]=size(Pest);[r4,~]=size(x_true);
if (r1~=1) | (r2~=1) | (r3~=1) | (r4~=1)
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
MCcolors.purple=[147,23,255]/255;
MCcolors.green=[0,160,0]/255;
MCcolors.red=[200,0,0]/255;
MCcolors.mag=[255,64,255]/255;
%
tbound=[t fliplr(t)];
bd=squeeze(Pest);
[len,wid]=size(bd);
if ((len>1) & (wid==1)) 
    bd=bd';
end
%
if strcmp(plot_type,'state')
    %bounds are around state estimate
    bound=[xest+2*sqrt(bd) fliplr(xest-2*sqrt(bd))];
    %
    plot(t,xest,'color',MCcolors.blue);
    hold on;
    patch(tbound,bound,'b','facecolor',MCcolors.blue,'edgecolor',MCcolors.blue,'linewidth',1,'FaceAlpha',.1,'EdgeAlpha',0.2);
    plot(t,x_true,'-.','Color',MCcolors.mag);
    if ~isempty(z)
        plot(t,z,'.','Color',MCcolors.green,'MarkerSize',12);
        legend('state estimate','estimate +/- 2{\sigma}','true state','measurement','Location','best');
    else
        legend('state estimate','estimate +/- 2\sigma','true state','Location','best');
    end
    hold off
    grid;
    xlabel(time_label);
    ylabel([['state estimate \it{x}' char(770)],'\it{(t)}'], 'Interpreter','tex', 'FontName','Helvetica')

end
%
if strcmp(plot_type,'error')
    err = xest - x_true;
    %bounds are around 0
    bound=[+2*sqrt(bd) fliplr(-2*sqrt(bd))];
    %
    plot(t,err,'-','Color',MCcolors.purple);
    hold on;
    patch(tbound,bound,'b','facecolor',MCcolors.blue,'edgecolor',MCcolors.blue,'linewidth',1,'FaceAlpha',.1,'EdgeAlpha',0.2);
    if ~isempty(z)
        plot(t,z-x_true,'.','Color',MCcolors.red,'MarkerSize',12);
        legend('error estimate','+/- 2\sigma bound','msmt noise','Location','best');
    else
        legend('error estimate','+/- 2\sigma bound','Location','best');
    end
    grid;
    xlabel(time_label);
    ylabel('error estimate {\it{e(t)}}');
end

end
