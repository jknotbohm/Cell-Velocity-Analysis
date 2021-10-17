function vel_autocorr_nogrid(fname)
%
% Autocorrelation of cell velocity data. The data do not have to be on a
% uniform grid and do not need to cover the entire field of view. This
% script requires a mat file containing cell trajectories.
%
% If running as a batch, uncomment the statement function at the top,
% comment the clear command, and comment the user input 'fname'
%
% Notes:
% - The correlation is performed on differences of displacements along the
%   cell trajectories. To compute velocity we'd divide by time, but doing
%   so has no effect on the autocorrelation, so it's not performed in this
%   script
% - This requires the functions autocorr_nogrid.m and
%   autocorr_angle_nogrid.m, available in the repository
%   'plotting_subfunctions'
% - The subfunction to correlate a unit vector can perform a polar
%   correlation (considering direction) or a nematic (considering only
%   orientation, with no distinction between forward and backward) This
%   script only keeps the polar correlation.
% - Take a look at the subfunctions to determine how the data are binned.
%   As of now, the bins are chosen such that they're centered at
%   spacing, 2*spacing, 3*spacing, etc. which is why variable r (in the
%   section for plotting) is defined as r = (0:1:Npts)*spacing.
% - This script can run on a time lapse data set
% - This code was originally written for a velocity autocorrelation for
%   data that's mapped to each cell (by computing cell trajectories). It
%   would also work for standard velocity data sets from DIC, though the
%   data would have to be pre-processed into the correct format. See notes
%   in the section LOAD DATA for info on formatting requirements.
% - There are two ways to run a correlation of velocity data--either the
%   vectors with magnitude can be correlated or the unit vectors (ie,
%   directions) can be correlated. The choice depends on your objective.
%   There's a user input to pick which option you're using.
% - There's an example of a way to plot the data included. You may want to
%   write a separate plotting script for your own purposes.
% - At the end is code to compute a correlation distance. The code computes
%   the correlation distance based on the mean autocorrelation over all
%   time points. If you were interested in changes in correlation distance
%   over time, you'd have to move the code for correlation distance into
%   the for loop.

%
% Written by Jacob Notbohm, University of Wisconsin-Madison 2020-2021
%


% clear;
close all;
clc;

%% --- USER INPUTS ---

% Name of mat file containing cell trajectories and tractions
% fname = 'cell_trajectories.mat';

% Number of points at which to compute the autocorrelation. This gives the
% length of the autocorrelation vector. Obviously, fewer points means the
% computation will run faster. But note that the relationship between
% computation and Npts is nonlinear--for larger distances (greater Npts),
% there are often more data points to correlate, causing more time than a
% correlation over a short distance.
Npts = 25;
% Spacing between the bins in the autocorrelation vector. If this is too
% small (with too small being smaller than the spatial resolution), the
% smallest distances will have nothing to correlate or very few points to
% correlate, which can lead to unexpected results, like 0s or negative
% numbers.
spacing = 20; % um
% Choose whether you're correting vectors with magnitudes or unit vectors.
corr_method = 2; % 1: vectors with magnitudes; 2: unit vectors


% Minimum value of correlation to show on plot axes. (The maximum value
% isn't a user input, because it's 1.)
corr_min = -0.2;
% Name to save plot. Set to empty array [] to suppress saving plot.
savename_plot = [];
% Name to save data. Set to empty array [] to suppress saving data.
savename_data = 'vel_autocorr_nogrid';

% After running the autocorrelations, the script computes a
% correlation distance, based on the distance at which a certain value of
% the correlation is reached. Typical values of corr_val range from 0.1 to
% 0.5. Before choosing the correct corr_val, look at the plots of the
% autocorrelation and verify that the spacing is chosen appropriately. Set
% to empty array [] to suppress finding the correlation distance.
corr_val = 0.3; % Unitless


%% --- LOAD DATA ---

% Load data
load(fname);
% Units: um
% Rows are the different cells
% Columns are different time points

%% --- RUN AUTOCORRELATION FOR EACH TIME POINT ---

% Number of time points
K = size(traj_x,2);

% Preallocate vector for all time points
AC_vel = zeros(K-1,Npts+1); % rows are the time points

for k=1:K-1
    % Get k-th cell positions
    traj_xk = traj_x(:,k);
    traj_yk = traj_y(:,k);
    % Get k-th cell velocity
    v_xk = traj_x(:,k+1)-traj_x(:,k);
    v_yk = traj_y(:,k+1)-traj_y(:,k);
    
    % Run autocorrelation
    AC_vel_k = autocorr_nogrid(traj_xk, traj_yk, v_xk, v_yk, Npts, spacing);
    % Run autocorrelation
    switch corr_method
        case 1
            % For vectors with magnitude
            AC_trac_k = autocorr_nogrid(traj_xk, traj_yk, v_xk, v_yk, Npts, spacing);
        case 2
            % For unit vectors. Note that the subfunction converts the vector to a unit vector.
            [AC_trac_k, ~] = autocorr_angle_nogrid(traj_xk, traj_yk, v_xk, v_yk, Npts, spacing);
        otherwise
            error('corr_method must be either 1 or 2; see description in USER INPUTS')
    end
    
    %     % Plot for checking/debugging. (Comment this when running all time points)
    %     figure;
    %     rr = (0:1:Npts)*spacing + spacing/2;
    %     plot(rr, AC_vel_k);
    %     xlabel('r (\mum)')
    %     ylabel('C(r)');
    %     set(gca,'box','off');
    
    
    % Add to array containing all data
    AC_vel(k,:) = AC_vel_k;
    
end



%% --- PLOT ALL DATA ---

hf = make_fig([0.2 0.2 1.2 1.2]);
r = (0:1:Npts)*spacing;
hold on
cmap = colormap('winter');
for k=1:size(AC_vel,1)
    ck = cmap(round(k/size(AC_vel,1)*size(cmap,1)) ,:);
    plot(r, AC_vel(k,:), 'color', ck);
end
plot(r, nanmean(AC_vel,1), 'k', 'linewidth', 2)
ylim([corr_min, 1])

xlabel('r (\mum)')
ylabel('C(r)');
set(gca,'box','off');

% Save plot
if ~isempty(savename_plot)
    set('paperpositionmode','auto');
    print(hf,'-dpng','-r300',savename_plot);
    % Option to save as eps
    % print(hf,'-depsc',savename_plot);
end

%% --- SAVE DATA ---
if ~isempty(savename_data)
    save(savename_data,'AC_vel','r','Npts','spacing')
end

%% --- GET CORRELATION DISTANCE ON MEAN OF TRACTIONS ---

AC_vel_mean = nanmean(AC_vel,1);

% Get correlation distance
if ~isempty(corr_val)
    % Get correlation size r_{corr_val}. This uses a cubic interpolation of all
    % data points surrounding the point corr_val
    idx = find(AC_vel_mean<corr_val, 1, 'first');
    
    if ~isnan(idx) & idx<=(length(AC_vel_mean)-1)
        if idx>2 && all(~isnan(AC_vel_mean( (-2:1)+idx )))
            corr_dist = interp1(AC_vel_mean( (-2:1)+idx ), r( (-2:1)+idx ), corr_val, 'pchip');
        elseif idx==2 && all(~isnan(AC_vel_mean( (-1:1)+idx )))
            corr_dist = interp1(AC_vel_mean( (-1:1)+idx ), r( (-1:1)+idx ), corr_val, 'pchip');
        else
            corr_dist = nan;
        end
    else
        corr_dist=nan;
    end
    
    % Save correlation distance
    writematrix(corr_dist,'vel_corr_dist.txt');
end
