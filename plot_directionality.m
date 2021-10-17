function plot_directionality(trajname)
% 
% This script gets directionality in a wound healing assay, which is
% defined as the cosine of the angle difference between the angle of 
% migration and the angle made by the unit normal vector of the wound edge.
%
% If running as a batch, uncomment the statement function at the top and
% comment the clear command
% 
% Notes
%   - The angle of the wound is assumed to be vertical, and the wound is 
%     assumed to be centered in the images, which may not match all
%     experimental images
%
% Written by Jacob Notbohm, University of Wisconsin-Madison, 2021
% 

% clear;
close all;
clc;

%% --- USER INPUTS ---

% % Name of data containing cell trajectories
% trajname = 'cell_trajectories.mat';

% Choose time points to use in the analysis. Select timepoints as a
% fraction as the total number of time points available, where 0
% corresponds to the first time point and 1 corresponds to the last time
% point.
nstart = 0;
nend = 1;

% Max directionality ratio - for color plots
directionality_max = 1;
% Set axis limits. To set the limits automatically, set to empty array [].
% Format [xmin xmax ymin ymax]. Units: um
xy_lim = []; 

% Name to save plot 
savename_plot = 'Directionality';
% Name to save data 
savename_data = 'Directionality_Data.mat';


%% --- LOAD DATA ---

load(trajname); % Units: um
% Rows indicate different trajectories
% Columns correspond to different time points

% --- Get time points of interest ---
K = size(traj_x,2);
idx = round(nstart*K)+1 : round(nend*K);
traj_x = traj_x(:,idx);
traj_y = traj_y(:,idx);

%% --- COMPUTE DIRECTIONALITY ---

% For each trajectory, get unit vectors corresponding to cell velocities at
% each time point
dx = traj_x(:,2:end) - traj_x(:,1:end-1);
dy = traj_y(:,2:end) - traj_y(:,1:end-1);
mag = sqrt(dx.^2+dy.^2);
dx_u = dx./mag;
% dy_u = dy./mag; % not needed

% Get directionality from dot product with unit vector in the horizontal 
% direction. The dot product is simply equal to dx_u.
directionality = dx_u;

% Determine if starting point is on left or right side of the image
xmin = min(traj_x(:,1));
xmax = max(traj_x(:,1));
% Indices correspond to cells on the right side of the image
I_right = traj_x(:,1) > (xmin+xmax)/2;
% If on the right side of the image, the dot product should be with unit
% vector in the -x direction, so switch the sign.
directionality(I_right,:) = -directionality(I_right,:);





%% --- PLOT RESULTS ---

hf1 = make_fig([0.5 1 2 .6]);

% Plot mean for each trajectory, given as a dot at starting
% location of each trajectory
subplot(1,3,1)
plot_colordot(traj_x(:,1),traj_y(:,1),nanmean(directionality,2),0,directionality_max,'jet')
caxis([0, directionality_max]);colorbar;
xlabel('(\mum)');
ylabel('(\mum)');
set(gca,'color','k');
title('Directionality');
% Axis limits
if ~isempty(xy_lim)
    axis(xy_lim);
else
    % Otherwise, get the automatically
    x1 = min(traj_x(:));
    x2 = max(traj_x(:));
    y1 = min(traj_y(:));
    y2 = max(traj_y(:));
    axis([x1 x2 y1 y2])
end
daspect([1 1 1]);

% Plot mean for each trajectory, given as a dot at final
% location of each trajectory
subplot(1,3,2)
plot_colordot(traj_x(:,end),traj_y(:,end),nanmean(directionality,2),0,directionality_max,'jet')
caxis([0, directionality_max]);colorbar;
xlabel('(\mum)');
ylabel('(\mum)');
set(gca,'color','k');
title('Directionality');
% Axis limits
if ~isempty(xy_lim)
    axis(xy_lim);
else
    % Otherwise, get the automatically
    x1 = min(traj_x(:));
    x2 = max(traj_x(:));
    y1 = min(traj_y(:));
    y2 = max(traj_y(:));
    axis([x1 x2 y1 y2])
end
daspect([1 1 1]);

% Histogram of all trajectories and all time points
subplot(1,3,3)

edges = linspace(0, directionality_max, 20);
histogram(directionality(~isnan(directionality)),edges,'Normalization','pdf');
set(gca,'box','off');
xlabel('Directionality')
ylabel('Approximate pdf')

% Save figure
print('-dpng','-r300',savename_plot);
% % Option to save as eps
% print('-depsc',savename_plot);
% Option to save in a different directory. Need to make sure that the
% selected path exists. Uncomment the next two lines and adjust the path as
% desired
% curdir = pwd;
% print('-dpng','-r300',['../Directionality/pos',curdir(end-1:end)]);

%% --- SAVE MSD EXPONENT AND PATH RATIO DATA ---

save(savename_data,'directionality');



