function plot_PathRatio(trajname)
% 
% Plot PathRatio, defined as end-to-end distance of the trajectory divided 
% by total distance. Cells that are more persistent have a larger path
% ratio.
% 
% Note that this definition differs from prior versions of the code, which 
% defined the path ratio to be the inverse of the currenent definition. 
% Before running this script, run compute_cell_trajectories.m.
%
% If running as a batch, uncomment the statement function at the top, 
% comment the clear command, and comment the user input 'trajname'
% 
% Notes
%   - The time between images is loaded by a separate file called time
%     increment.
%
% Written by Jacob Notbohm, University of Wisconsin-Madison, 2014-2021
% 

% clear;
close all;
clc;

%% --- USER INPUTS ---

% Name of data containing cell trajectories
% trajname = 'cell_trajectories.mat';

% Choose time points to use in the analysis. Select timepoints as a
% fraction as the total number of time points available, where 0
% corresponds to the first time point and 1 corresponds to the last time
% point.
nstart = 0;
nend = 1;

% Max path length ratio - for color plots
pratio_max = 1;
% Set axis limits. To set the limits automatically, set to empty array [].
% Format [xmin xmax ymin ymax]. Units: um
xy_lim = []; 

% Name to save plot 
savename_plot = 'PathRatio';
% Name to save data 
savename_data = 'PathRatio_Data.mat';


%% --- LOAD DATA ---

% Get time between images
% % copyfile('../TimeIncrement.txt','TimeIncrement.txt');
fid = fopen('TimeIncrement.txt');
txtcell = cell2mat(textscan(fid,'%f %*[^\n]')); % '%*[^\n]' skips the remainder of each line
time_increment = txtcell(1); % min
fclose(fid);

load(trajname); % Units: um


%% --- COMPUTE PATH RATIO: PATH LENGTH / END-TO-END DISTANCE ---

% --- Get time points of interest ---
K = size(traj_x,2);
idx = round(nstart*K)+1 : round(nend*K);
traj_x = traj_x(:,idx);
traj_y = traj_y(:,idx);

% Only use trajectories that don't have a nan component
idx_x = ~any(isnan(traj_x),2); % index of rows with no nans
idx_y = ~any(isnan(traj_y),2);
idx = idx_x & idx_y;
traj_x = traj_x(idx,:);
traj_y = traj_y(idx,:);

% End to end distance
D_endtoend = sqrt( (traj_x(:,end)-traj_x(:,1)).^2 + (traj_y(:,end)-traj_y(:,1)).^2 );
% Incremental pathlength
D_incremental_x = traj_x(:,2:end)-traj_x(:,1:end-1);
D_incremental_y = traj_y(:,2:end)-traj_y(:,1:end-1);
D_incremental = sqrt( D_incremental_x.^2 + D_incremental_y.^2 );
% Cumulative path length
D_path = sum(D_incremental,2);

% Path ratio
path_ratio_all =  D_endtoend ./ D_path;
% Mean path ratio
path_ratio = mean(path_ratio_all);


%% --- PLOT RESULTS ---

hf1 = make_fig([0.5 1 1.2 .6]);

% Second plot is pathlength for each trajectory, given as a dot at starting
% location of each trajectory
subplot(1,2,1)
plot_colordot(traj_x(:,1),traj_y(:,1),path_ratio_all,0,pratio_max,'jet')
caxis([0, pratio_max]);colorbar;
xlabel('(\mum)');
ylabel('(\mum)');
set(gca,'color','k');
title('Path length ratio');
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

% Third plot is pathlength for each trajectory, given as a dot at ending
% locating of each trajectory
subplot(1,2,2)
plot_colordot(traj_x(:,K),traj_y(:,K),path_ratio_all,0,pratio_max,'jet')
caxis([0, pratio_max]);colorbar;
xlabel('(\mum)');
ylabel('(\mum)');
set(gca,'color','k');
title('Path length ratio');
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


% Save figure
print('-dpng','-r300',savename_plot);
% % Option to save as eps
% print('-depsc',savename_plot);
% Option to save in a different directory. Need to make sure that the
% selected path exists. Uncomment the next two lines and adjust the path as
% desired
% curdir = pwd;
% print('-dpng','-r300',['../MSD/pos',curdir(end-1:end)]);

%% --- SAVE MSD EXPONENT AND PATH RATIO DATA ---

save(savename_data,'path_ratio_all');
% Outputs:
%   path_ratio_all: path ratio for each trajectory


