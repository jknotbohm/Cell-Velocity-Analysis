function plot_MSD(trajname)
% 
% Mean square displacement following trajectories from DIC data. Before
% running this script, run compute_cell_trajectories.m.
%
% If running as a batch, uncomment the statement function at the top, 
% comment the clear command, and comment the user input 'trajname'
% 
% Notes
% - The time between images is loaded by a separate file called time
%   increment.
% - After the MSD is computed, there are some commented lines for
%   subanalyses to compute MSD a certain distance from the center. These
%   lines of code assume that the origin (0,0) is the center of the cell
%   layer, which will be the case if the script compute_cell_trajectories
%   shifts the data so that the origin is at the center--this is typically
%   done for experiments with a cell island. These lines of code are
%   representative of subanalyses that can be performed. Another option is
%   to compute MSD for a wound healing experiments, where cells are grouped
%   by distance from the free edge--such an ananlysis isn't yet coded, but
%   it could be done using a domain image to determine the free edge.
%
% Written by Jacob Notbohm, University of Wisconsin-Madison, 2014-2021
% 

% clear;
close all;
clc;

%% --- USER INPUTS ---

% % Name of data containing cell trajectories
trajname = 'cell_trajectories.mat';

% Choose time points to use in the analysis. Select timepoints as a
% fraction as the total number of time points available, where 0
% corresponds to the first time point and 1 corresponds to the last time
% point.
nstart = 0;
nend = 1;

% Name to save plot 
savename_plot = 'MSD_plot';
% Name to save data 
savename_data = 'MSD_Data.mat';


%% --- LOAD DATA ---

% Get time between images
% % copyfile('../TimeIncrement.txt','TimeIncrement.txt');
fid = fopen('TimeIncrement.txt');
txtcell = cell2mat(textscan(fid,'%f %*[^\n]')); % '%*[^\n]' skips the remainder of each line
time_increment = txtcell(1); % min
fclose(fid);

load(trajname); % Units: um

%% --- COMPUTE MSD ---

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

% Number of trajectories
N = size(traj_x,1);
% Number of time points
K = size(traj_x,2);

% Compute MSD for each time shift, dt
MSD = zeros(N,K-1)*nan; % Preallocate
for dt = 1:K-1 % Units: time points
    for n=1:N
        % Positions of trajectory n at reference time t0
        x1 = traj_x(n,1:K-dt);
        y1 = traj_y(n,1:K-dt);
        % Positions of trajectory n at time t0+dt
        x2 = traj_x(n,(1:K-dt)+dt);
        y2 = traj_y(n,(1:K-dt)+dt);
        % Mean square displacement for n-th trajectory and time shift dt
        MSD(n,dt) = mean((x1-x2).^2 + (y1-y2).^2);
    end
end

% Mean MSD for all trajectories
MSD_mean = mean(MSD,1);

% % Mean MSD for trajectories starting within 250 um of center
% x1 = traj_x(:,1);
% y1 = traj_y(:,1);
% dist = sqrt(x1.^2+y1.^2);
% idx = dist<250;
% MSD_mean_bulk = mean(MSD(idx,:),1);
% 
% % Mean MSD for trajectories starting outside 250 um of center
% idx = dist>250;
% MSD_mean_edge = mean(MSD(idx,:),1);
% % Create vector of time shifts dt
% dt = (1:K-1)*time_increment; % units: min

% % Mean MSD for trajectories within 1600 um from center
% x1 = traj_x(:,1);
% y1 = traj_y(:,1);
% dist = sqrt(x1.^2+y1.^2);
% idx = dist<1600;
% MSD_mean = mean(MSD(idx,:),1);


%% --- PLOT ---

dt = (1:K-1)*time_increment; % units: min

hf1 = make_fig([0.5 1 .6 .6]);

hold on
plot(dt,MSD_mean,'b','linewidth',2);
% Also plot lines with slope of 1 and 2 in log scale
y_slope1 = dt * MSD_mean(1)/dt(1); % Add multiplier to make initial points match
plot(dt,y_slope1,'color',[0.7 0.7 0.7]);
y_slope2 = dt.^2 * MSD_mean(1)/dt(1)^2;
plot(dt,y_slope2,'color',[0.7 0.7 0.7]);
set(gca,'xscale','log','yscale','log');
xlabel('\Deltat (min)');
ylabel('MSD(\Deltat) (\mum^2)');

% Fit to power law
p = polyfit(log10(dt),log10(MSD_mean),1);
n = p(1);

% Write power on plot
annotation(hf1,'textbox',[0.15 0.9 0.6 0.03],...
    'String',{['Slope = ',num2str(n,'%3.2f')]},'FitBoxToText','off','fontsize',11,'linestyle','none');
% % Write directory on plot
% annotation(hf1,'textbox',[0.003 0.021 0.98 0.03],...
%     'String',{pwd},'FitBoxToText','off','fontsize',8,'linestyle','none','interpreter','none');

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

save(savename_data,'n','dt','MSD_mean');
% Outputs:
%   n: MSD exponent
%   dt: MSD delta t (min)
%   MSD_mean: mean MSD of cell trajectories (um^2)

