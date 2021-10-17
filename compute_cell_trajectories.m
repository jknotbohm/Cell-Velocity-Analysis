function compute_cell_trajectories
% compute_cell_trajectories
%
% Compute pathlines for cell trajectories based on cell displacements
% calculated from image correlation.
%
% Algorithm adapted from a code written by Chan Young Park, Harvard School
% of Public Health, 2012
% 
% If running as a batch, uncomment the statement function at the top and
% comment the clear command
% 
% Several important notes:
% - A domain file is needed for at least the first time point. It's possible
%   that results will be more accurate if you have a domain file for all
%   time points. Get the domain file by running find_boundaries.m
% - The file ExperimentalSettings.txt, which is used in running traction
%   analysis is required. This file is used to get the pixel size.
% - There's a section that performs drift correction by subtracting the
%   median, which may or may not work for your application. Take a look at 
%   the comments, and test whether this improves your data. By default, the
%   median subtraction is commented out.
% 
% Written by Jacob Notbohm, University of Wisconsin-Madison, 2014-2021

% clear;
close all;
clc;


%% --- USER INPUTS ---
% Name of displacment data to load
DICname = 'DIC_results_w0=64_inc.mat';
% Name of domain. This is where cells are located. The domain comes from 
% running find_boundary.m. Set to [] if no domain
domainname = 'domain01.tif';

% Factor by which to downsample data. (Need to be careful here -- final
% number of grid points should match number of cells.) Must be positive
% integer.
fd = 2;

% State whether geometry is an island. Set this to 1 if yes. For island
% geometry, code will set origin to be the center of the island
isisland = 0;
% Choose starting and ending time points for analysis. Enter
% empty array [] to begin at first time point and/or end at last time point
% tstart = 1;
% tend = [];
tstart_end = load('time_points_start_end.txt');
tstart = tstart_end(1);
if isnan(tstart_end(2))
    tend = [];
else
    tend = tstart_end(2);
end

% Threshold to reject spurious displacement data. Any incremental 
% displacement data greater than this value is set to nan. Use a large 
% value to avoid thresholding.
% Note that units here are um, whereas in some other code units are um/min
thr = 15; % um

% Name to save data
savename = 'cell_trajectories_tstart_end.mat';


%% --- GET TRAJECTORIES ---

% Get pixel size from experimental settings file
fid = fopen('ExperimentalSettings.txt');
txtcell = cell2mat(textscan(fid,'%f %*[^\n]')); % '%*[^\n]' skips the remainder of each line
pix_size = txtcell(1)*1e6; % Pixel size, um

% Load data
load(DICname);
% Rename variables associated with cell displacements so they aren't
% confused with substrate displacements
u_cell=u; v_cell=v; x_cell=x; y_cell=y; w0_cell=w0; d0_cell=d0;
% Convert from pix to um
x_cell=x_cell*pix_size;     y_cell=y_cell*pix_size;
u_cell=u_cell*pix_size;     v_cell=v_cell*pix_size;

% Check if domain is entire image
if isempty(domainname)
    info = imfinfo(cellname);
    domain1 = ones(info(1).Height,info(1).Width);
else % Otherwise load domain
    domain1 = imread(domainname,1);
    domain1 = double(domain1); % Convert to double precision
    domain1 = domain1/max(domain1(:)); % Set max value to 1
end

% Downsample domain
% x and y grid points don't necessarily start at 1.
% First crop off edges so that domain matches start and end points of x
% and y.
domain1 = domain1(min(y(:)):max(y(:)), min(x(:)):max(x(:))); % Match starting row and col to y and x
domain1 = downsample(domain1,d0*fd); % downsample number of rows
domain1 = downsample(domain1',d0*fd)'; % downsample number of cols

% If geometry is island, shift origin so it's at center of domain
if isisland == 1
    % Get center from center of first domain image
    % Centroid coordinates
    xc = sum(x(:).*domain1(:)) / sum(domain1(:)); % Units: pix (same units as x and y)
    yc = sum(y(:).*domain1(:)) / sum(domain1(:));
    % Convert to um
    xc = xc*pix_size;
    yc = yc*pix_size;
    
    % Shift data so origin is at center of domain
    x_cell = x_cell - xc; % Units: um
    y_cell = y_cell - yc;
end

% Get initial grid points for trajectory computation
x_cell2=downsample(x_cell,fd); x_cell2=downsample(x_cell2',fd)';
y_cell2=downsample(y_cell,fd); y_cell2=downsample(y_cell2',fd)';



% Get domain for first time point
% For an island, find the largest connected region
if isisland == 1
    BNDRY = bwboundaries(domain1); % This should be 1 cell with boundary coordinates
    LB = cellfun(@length,BNDRY);
    [~, idx] = max(LB);
    BNDRY = BNDRY{idx};
    x_bndry = BNDRY(:,2); % x coordinates are the columns
    y_bndry = BNDRY(:,1); % y coordinates are the rows
    % Convert from d0-spaced pix to um
    x_bndry=x_bndry*pix_size*d0; y_bndry=y_bndry*pix_size*d0;
    % Shift so origin is at zero
    x_bndry = x_bndry-xc; % Units: um
    y_bndry = y_bndry-yc;
    
    % Get indices of points inside region given by (x_bndry,y_bndry)
    IDX = inpolygon(x_cell2,y_cell2,x_bndry,y_bndry);
    
% For edge geometry, just use the given domain    
else
    IDX = logical(domain1);
end


% Determine starting time point for analysis. If first images are all
% zeros, then starting timepoint should be 2.
if isempty(tstart)
    u1 = u_cell(:,:,1);
    v1 = v_cell(:,:,1);
    if sum(u1(:))==0 && sum(v1(:))==0
        tstart=2;
    else
        tstart=1;
    end
end

% Determine final time point for analysis
if isempty(tend)
    tend = size(u_cell,3)-1;
end
% Make sure tend is smaller than the 3rd dimension of u_cell. If it's
% larger by 1, it may be because tend correspond to the last time point,
% but there's one fewer correlation than time points, in which case we can
% simply subtract 1 from tend. Otherwise, there's a problem with the value
% of tend input
if tend==size(u_cell,3)+1
    tend = size(u_cell,3);
elseif tend>size(u_cell,3)
    error('tend is larger than the number of time points in the DIC data');
end

% Number of trajectories is given by number of nonzero elements in IDX
num_traj = sum(double(IDX(:)));

% Preallocate arrays for trajectories. These arrays start at each 
% downsampled  grid point within the domain. Rows are different traj and 
% columns are different times.
traj_x = nan*zeros(num_traj,tend-tstart+1);
traj_y = nan*zeros(num_traj,tend-tstart+1);

% First set of trajectories are given by points inside domain (ie, where IDX==1)
traj_x(:,1) = x_cell2(IDX);
traj_y(:,1) = y_cell2(IDX);

% Scan through remaining time points adding to the path coordinates
for k= tstart:tend
    if ~isempty(domainname)
        
        % Determine whether domain is for only first time point or for all
        % time points
        info = imfinfo(domainname);
        
        % If there's a domain for each time point, then get k-th domain
        if length(info)>1
            % There's an ambiguity here, because k is correlation number, 
            % and the correlation is between the k-th and (k+1)-th images 
            % if tstart==1 or the (k-1)-th and k-th images if tstart==2.
            % Here, I'm choosing the starting image as k-tstart+1
            domain = imread(domainname, k-tstart+1);
            
            domain = double(domain);
            domain = domain/max(domain(:));
            domain = logical(domain);
            % Downsample domain
            % x and y grid points start at w0/2 and end w0/2 before the image ends.
            % First crop off edges so that domain matches start and end points of x
            % and y.
            domain = domain(min(y(:)):max(y(:)), min(x(:)):max(x(:)));
            domain = downsample(domain,d0); % downsample number of rows
            domain = downsample(domain',d0)'; % downsample number of cols
        else
            % Otherwise, keep all data in the image
            domain = true(size(x));
        end
    else
        % If no domain name, keep all data in the image
        domain = true(size(x));
    end
    
    % --- k-th cell displacements ---
    u_cell_k = u_cell(:,:,k); % units: um
    v_cell_k = v_cell(:,:,k);
    
    % --- Remove displacements that are too large ---
    % Also include a threshold based on size of image
    xymax = sqrt(max(x_cell(:))*max(y_cell(:))); % Typical image size. Units: um
    % Threshold for velocity data
    idx = (abs(u_cell_k)>thr) | (abs(v_cell_k)>thr) | abs(u_cell_k)>0.1*xymax | abs(v_cell_k)>0.1*xymax;
    u_cell_k(idx) = nan;
    v_cell_k(idx) = nan;
    u_cell_k = inpaint_nans(u_cell_k);
    v_cell_k = inpaint_nans(v_cell_k);
    
    % --- Correct for drift. Choose one of various options ---
    % Be carful -- These may or may not work for you.
    % Ideally, you'd determine drift by measuring a motion of a reference
    % point in the dish. With only the cell images, such a reference point
    % doesn't exist. The code below subtracts the median velocity, which
    % may or may not be a good idea. It's commented out by default. Perhaps
    % a better option would be to use images of fluorescent particles in a
    % TFM experiment to determine drift. I did this once in the past, and
    % didn't seem to get improved results. It's worth trying again.
    
%     % Correct for drift by subtracting off mean displacement of the slowest
%     % third of the cells
%     idx = ~isnan(u_cell_k) & ~isnan(v_cell_k);
%     u_tmp = u_cell_k(idx);
%     v_tmp = v_cell_k(idx);
%     umag_tmp = sqrt(u_tmp.^2+v_tmp.^2);
%     umag30 = prctile(umag_tmp,30);
%     idx = umag_tmp < umag30;
%     uave = mean( u_tmp(idx) );
%     vave = mean( v_tmp(idx) );
    
%     % Correct for drift by subtracting off medians
%     uave = nanmedian(u_cell_k(:));
%     vave = nanmedian(v_cell_k(:));

    % No drift correction [default]
    uave=0; vave=0;
    
    u_cell_k = u_cell_k-uave;
    v_cell_k = v_cell_k-vave;
    
    % Set data outside domain to nan
    u_cell_k(~domain)=nan;
    v_cell_k(~domain)=nan;
    
    % Interpolate k-th displacements to gridpoints of previous timepoint
    displ_x = griddata(x_cell,y_cell,u_cell_k,traj_x(:,k-tstart+1),traj_y(:,k-tstart+1));
    displ_y = griddata(x_cell,y_cell,v_cell_k,traj_x(:,k-tstart+1),traj_y(:,k-tstart+1));
        
    % Add to trajectory arrays
    traj_x(:,k-tstart+2) = traj_x(:,k-tstart+1) + displ_x; % Units: um
    traj_y(:,k-tstart+2) = traj_y(:,k-tstart+1) + displ_y;
end


%% --- IF THERE'S A TIME INCREMENT FILE, THEN ADD VELOCITY DATA ---

% Get time between images
if isfile('TimeIncrement.txt')
    fid = fopen('TimeIncrement.txt');
    txtcell = cell2mat(textscan(fid,'%f %*[^\n]')); % '%*[^\n]' skips the remainder of each line
    time_increment = txtcell(1); % min
    fclose(fid);
    
    dtx = traj_x(:,2:end)-traj_x(:,1:end-1);
    dty = traj_y(:,2:end)-traj_y(:,1:end-1);
    speed = sqrt(dtx.^2 + dty.^2)/time_increment; % um/min
    % Get trajectories for which a time point is nan and remove
    I = any(isnan(speed),2);
    speed2 = speed(~I,:);
    ave_speed = nanmean(speed2(:));
end


%% --- SAVE DATA ---
if isfile('TimeIncrement.txt')
    save(savename,'traj_x','traj_y','ave_speed'); % Units: um; um/min
else
    save(savename,'traj_x','traj_y'); % Units: um
end




