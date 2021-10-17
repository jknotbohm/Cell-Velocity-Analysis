function plot_D2min(trajname)
% 
% Compute D2min, where D2min is defined by Lee et al New J. Physics 15
% (2013) 025036 as
% 
% D2min(i) = min_E{ 1/N*sum_j[(Delta d_ij - E*d_ij)^2] }
% 
% where i indexes a cell and j indexes the nearest N neighbors of cell i.
% Often we pick cells indexed by j to be within some distance of cell i,
% say 40 um. The sum is over the j cells. The minimization is over E, which
% is a strain tensor describing the deformation of cell i and the cells
% indexed by j. d_ij is defined as the relative position vector connecting
% cells i and j. Delta d_ij is a change in d_ij over time. Generally the
% time is equal to the time required for a cell to move 1 cell diameter.
% 
% Before running this script, run compute_cell_trajectories.m.
%
% If running as a batch, uncomment the statement function at the top, 
% comment the clear command, and comment the user input 'trajname'
% 
% Notes:
% - The time between images is loaded by a separate file called time
%   increment.
% - The equation is written as a sum over all cells within a certain
%   distance of the cell of interst. A similar equation is to take the sum
%   over the nearest N cells. We take this second approach, which is
%   preferred, because in some experimental geometries, the cells get 
%   farther away from each other over time.
% - The equation for the strain E comes from ref. 29 of the manuscript by
%   Lee et al. See the comments near the derivation for more info.
% 
% Written by Jacob Notbohm, University of Wisconsin-Madison, 2015-2021
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

% Specify number of nearest neighbors
n_nearestneighbors = 8; % A common number to use is 8
% Time shift used to compute Delta d_{ij} in equation above. Often this is
% equal to the time required to move 1 cell diameter.
deltaT = 10; % Units: time increments

% Limit for color plots
D2min_lim = 150; % Units: um^2

% For long timelapse data, we often don't need to plot results for every
% time point. Select how many time points are skipped between plots that
% are made. Set this parameter to 1 to plot all images and to a large
% number to plot only the first image.
im_inc = 20;

% Name of a subfolder to create to save plots in
dirhead = ['D2min_n',num2str(n_nearestneighbors)];
% Header of name to save plots. This can contain a directory listing
savenameheader = [dirhead,'/t_'];
% Set to [] to make figures created in the for loop visible
invisible = [];
% Set axis limits. To set the limits automatically, set to empty array [].
% Format [xmin xmax ymin ymax]. Units: um
xy_lim = []; 


%% --- LOAD DATA ---

% Make diretory to save data
if exist(dirhead,'dir')==7
    rmdir(dirhead,'s');
end
mkdir(dirhead);

% Get time between images
% % copyfile('../TimeIncrement.txt','TimeIncrement.txt');
fid = fopen('TimeIncrement.txt');
txtcell = cell2mat(textscan(fid,'%f %*[^\n]')); % '%*[^\n]' skips the remainder of each line
time_increment = txtcell(1); % min
fclose(fid);

load(trajname); % Units: um
% Note: These trajectories are often computed from DIC data. Technically,
% we're supposed to compute D2min from directly tracking cells and not
% correlations of subsets. When possible, calculate D2min from trajectories
% based on tracking nuclei


%% --- COMPUTE D2MIN AND MAKE PLOTS ---

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

% Number of trajectories, N, and number of time points, K
[N, K] = size(traj_x);

% Preallocate matrix of values
D2min = zeros(N,K-deltaT)*nan;
for k=1:K-deltaT-1
    % Trajectory positions at time point k
    traj_x_k = traj_x(:,k);
    traj_y_k = traj_y(:,k);
    % Trajectory positions at time point k+. Note that
    % D2min is often computed over a period of time equal to the time
    % required for a cell to move one cell diameter. Call this time deltaT,
    % where deltaT is a specific number of time steps
    traj_x_kplus = traj_x(:,k+deltaT);
    traj_y_kplus = traj_y(:,k+deltaT);
    
    % Compute D2min for each trajectory
    for i=1:N % i indexes over cell trajectories
        % Coordinates of cell i
        xi = traj_x_k(i); % Units: um
        yi = traj_y_k(i);
        
        xi_plus = traj_x_kplus(i);
        yi_plus = traj_y_kplus(i);
        
        % Get cells within a radius D2min_rad of (xi,yi)
        dist = sqrt((traj_x_k-xi).^2 + (traj_y_k-yi).^2); % Dist is a vector
        %         idx = dist<D2min_rad & dist>0;
        
        % Instead of using a distance, specify number of nearest neighbors
        [dist_sort, idx] = sort(dist,'ascend'); % Sort distances
        idx = idx(1:n_nearestneighbors); % Get n nearest neighbors
        
        xj = traj_x_k(idx); % Units: um
        yj = traj_y_k(idx);
        
        xj_plus = traj_x_kplus(idx);
        yj_plus = traj_y_kplus(idx);
        
        % Compute Delta d_ij
        Delta_dx = (xj_plus-xi_plus) - (xj-xi);
        Delta_dy = (yj_plus-yi_plus) - (yj-yi);
        
        % Compute strain tensor E that best matches the motion of cell i
        % and the cells indexed by j. I compute E from Eqs. 2.12-2.14 of
        % Falk & Langer, PRE 57, 6, 1998 (ref. 29 of Lee et al). 
        % The notation below is slightly different from falk & Langer, but
        % matches more closely the notation of Lee et al
        X(1,1) = sum((xj_plus-xi_plus).*(xj-xi)); % Sum is over j, the nearest neighbors
        X(1,2) = sum((xj_plus-xi_plus).*(yj-yi));
        X(2,1) = sum((yj_plus-yi_plus).*(xj-xi));
        X(2,2) = sum((yj_plus-yi_plus).*(yj-yi));
        
        Y(1,1) = sum((xj-xi).*(xj-xi));
        Y(1,2) = sum((xj-xi).*(yj-yi));
        Y(2,1) = sum((yj-yi).*(xj-xi));
        Y(2,2) = sum((yj-yi).*(yj-yi));
        
        % E = X*inv(Y) - eye(2); % * should be matrix multiplication; eye(2) is 2x2 identity matrix
        E = X/Y - eye(2); % Better implementation of the matrix equation on the line above
        
        % There may seem to be a diference between Lee et al and Falk &
        % Langer. Falk & Langer use (E+eye(2)) rather than just E in the
        % equation. But note that the first term in their equation is
        % different, and doesn't include (xj-xi) or (yj-yi); including the
        % eye(2) makes the equations used in both manuscripts the same.
        
        % Compute "affine displacement" by computing strain E times d_ij
        Edx = E(1,1)*(xj-xi) + E(1,2)*(yj-yi);
        Edy = E(2,1)*(xj-xi) + E(2,2)*(yj-yi);
        
        % Take the mean over all indices j to get D2min(i)
        D2min(i,k) = mean( (Delta_dx-Edx).^2 + (Delta_dy-Edy).^2 );
        
        % Another way of writing. These look to be the same.
%         D2min(i,k) = mean( ( (xj_plus-xi_plus) - ( (E(1,1)+1)*(xj-xi)+E(1,2)*(yj-yi) ) ).^2 + ...
%             ( (yj_plus-yi_plus) - ( E(2,1)*(xj-xi)+(E(2,2)+1)*(yj-yi) ) ).^2 );
        
        
    end
	
    % Make colored plot
    if k/im_inc == round(k/im_inc)
        hf1 = make_fig([0.1 0.5 0.8 0.8]);
        if invisible
            set(hf1,'visible','off');
        end
        plot_colordot(traj_x(:,k),traj_y(:,k),D2min(:,k),0,D2min_lim,'hot')
        caxis([0 D2min_lim]);colorbar;
        xlabel('\mum');
        ylabel('\mum');
        set(gca,'color','k');
        title('D^2_{min}');
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
        
        % Write directory on plot
        annotation(hf1,'textbox',[0.003 0.025 0.98 0.03],...
            'String',{pwd},'FitBoxToText','off','fontsize',8,'linestyle','none','interpreter','none');
        
        % Save
        set(hf1,'PaperPositionMode','auto','InvertHardCopy','off');
        print('-dpng','-r150',[savenameheader,num2str(k,'%3.3d')]);
        close(hf1);
    end
	
end

%% --- PLOT ALL DATA ---

hf = make_fig([0.8 0.5 1.4 1.4]);

% --- Plot D2min for individual cells ---

% For each cell trajectory, plot D2_min over time. This is used to identify
% if some cells consistently exhibit higher values of D2min, which might be
% typical of a leader cell.
% Get time variables
time = (0:size(D2min,2)-1)*time_increment;
% Downsample D2min for plotting
num_traj = size(D2min,1);
D2min_plot = D2min(1:20:num_traj,:);

subplot(2,2,1)
plot(time',D2min_plot');
ylim([0 2*D2min_lim])
xlabel('Time (min)');
ylabel('D^2_{min} (\mum^2)');
set(gca,'box','off');
title('D^2_{min} for selected trajectories')

% Histogram
subplot(2,2,2)
edges = linspace(0,2*D2min_lim);
[f, ~] = histcounts(D2min(:),edges);
bins = edges(1:end-1) + 0.5*(edges(2)-edges(1));
f = f/length(D2min(:)) / (bins(2)-bins(1)); % Change to frequency
plot(bins,f*100);
xlabel('D^2_{min} (\mum^2)')
ylabel('%')
set(gca,'box','off');
title('Histogram')

% --- Plot mean ---

% Compute mean of D2min over time
D2min_mean = nanmean(D2min,2);

% D2min plotted on cell positions at first time point
subplot(2,2,3)
plot_colordot(traj_x(:,1),traj_y(:,1),D2min_mean,0,D2min_lim,'hot')
caxis([0 D2min_lim]);colorbar;
xlabel('\mum');
ylabel('\mum');
set(gca,'color','k');
title('D^2_{min}');
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

% D2min plotted on cell positions at first last point
subplot(2,2,4)
plot_colordot(traj_x(:,K-deltaT),traj_y(:,K-deltaT),D2min_mean,0,D2min_lim,'hot')
caxis([0 D2min_lim]);colorbar;
xlabel('\mum');
ylabel('\mum');
set(gca,'color','k');
title('D^2_{min}');
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



set(hf,'Paperpositionmode','auto','InvertHardCopy','off');
print(hf,'-dpng','-r300',['D2min_mean_n',num2str(n_nearestneighbors)']);


%% --- SAVE DATA ---

save(['D2min_t06_n',num2str(n_nearestneighbors),'.mat'],'D2min','D2min_mean');

