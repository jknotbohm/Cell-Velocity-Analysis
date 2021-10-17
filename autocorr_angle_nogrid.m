function [ACp, ACn] = autocorr_angle_nogrid(varargin)
%AUTOCORR_ANGLE_NOGRID Autocorrelation of angle (both polar and nematic)
%
% This function performs the correlation by binning data and computing the
% autocorrelation of all data within each bin.
%
% This function runs autocorrelation of the direction of a vector.
%
% INPUTS
% x, y      vectors identifying points
% Tx, Ty    vectors identifying variable to autocorrelate (these can be
%           unit vectors, but they don't have to be; the code converts them
%           to unit length)
% Npts      Number of points at which to compute the autocorrelation. This
%           gives the length of the autocorrelation vector
% spacing   Spacing between the bins in the autocorrelation vector
% CellIdx   (Optional parameter) Scalar corresponding to index of the cell
%           for which to compute the autocorrelation. If no value is input,
%           the autocorrelation is computed and averaged for all cells,
%           which is likely the more common usage of this function.
%
% OUTPUTS
% ACp       Polar autocorrelation
% ACn       Nematic atocorrelation
%
% Note: The output should start with a value of 1, ie, ACp(1)==1 and
% ACn(1)==1. Subsequent values should not be NaN; if they are NaN, then no
% points exist within the requested spacing--increase the spacing.
%
% Written by Jacob Notbohm, University of Wisconsin-Madison 2020-2021
%

% Check outputs
if nargout ~= 2
    error('Wrong number of function outputs')
end

% Parse inputs
x = varargin{1};
y = varargin{2};
Tx = varargin{3};
Ty = varargin{4};
Npts = varargin{5};
spacing = varargin{6};
if length(varargin)==6
    F = 1;
elseif length(varargin)==7
    CellIdx = varargin{7};
    F = 2;
else
    error('Wrong number of inputs')
end

% Check for nan values
switch F
    case 1 % When doing autocorr over all data, can simply remove nan values
        idx = ~isnan(Tx) & ~isnan(Ty);
        x=x(idx); y=y(idx);
        Tx=Tx(idx); Ty=Ty(idx);
    case 2
        % If there are nan values in case 2, then they must be removed
        % before calling this function, because removing nans will mess up
        % CellIdx
        if any(isnan(Tx)) || any(isnan(Ty))
            error('When specifying CellIdx, there can be no nan values in Tx and Ty.')
        end
end

% Conver to unit vectors
Tmag = sqrt(Tx.^2+Ty.^2);
Tx2 = Tx./Tmag;
Ty2 = Ty./Tmag;

% Get angle
theta = atan2(Ty2, Tx2);

% --- Identify edges to discretize data ---
% There are a couple methods to do this:
% Method 1 (not currently in use). Start at 0 and add an edge at every
% increment of spacing.
% edges = (0:1:Npts)*spacing; % um
% % D has many 0 values that we want to ignore. Accomplish this by
% % setting lower limit of edges to a value slightly greater than 0.
% edges(1) = 1e-6;
% Method 2. Start at spacing/2. I prefer this method, because it means that
% the centers of the bins are given by spacing, 2*spacing, 3*spacing, etc.
edges = (0:1:Npts)*spacing + spacing/2; % um


% --- Run autocorr for the two different cases ---

switch F
    
    case 1 % --- AUTOCORR FOR ALL CELLS AND TAKE AVERAGE ---
        
        % For each point, get distance between that point and all other points
        N = length(x);
        D = zeros(N,N);
        for n=1:N
            D(:,n) = sqrt( (x-x(n)).^2 + (y-y(n)).^2 );
        end
        % Since D is symmetric, only need upper triangular portion. This will
        % avoid double counting cells and will cut computation time in half
        D = triu(D);
        
        % Discretize data
        I = discretize(D, edges);
        % nan values occur for cells that are a distance beyond the max of edges
        
        % --- Perform correlation ---
        
        % For each bin, find all cells in that bin and perform the correlation
        ACp = zeros(1, length(edges)-1)*nan; % preallocate
        ACn = ACp;
        for n=1:length(edges)-1
            I_n = find(I==n);
            % Get the rows and columns of I corresponding to I_n. The rows and
            % columns are indices for each pair of cells that we want to
            % correlate
            [I1, I2] = ind2sub([N,N] ,I_n);
            
            
            % Make sure there is more than one pair of points to correlate
            if length(I_n) <= 1
                AC(n) = nan;
            else
                % Perform autocorrelation
                % numerator
                theta1 = theta(I1);
                theta2 = theta(I2);
                
                % Polar autocorrelation. This is mathematically equivalent
                % to an autocorrelation of the unit vectors. This can be shown
                % using the trig identity cos(a-b) = cos(a)cos(b) + sin(a)sin(b).
                ACp_n = sum(cos(theta1-theta2));
                % Nematic autocorrelation.
                ACn_n = sum(cos(2*(theta1-theta2)));
                % Normalize by sum of cos(theta1-theta1) for all points, ie, by the
                % number of points
                ACp(n) = ACp_n / numel(theta1);
                ACn(n) = ACn_n / numel(theta1);
                
            end
        end
        
        % Autocorrelation for shift of 0 is equal to 1
        ACp = [1, ACp];
        ACn = [1, ACn];
        
        %     % Plot for checking/debugging
        %     figure;
        %     rr = (0:Npts)*spacing + spacing/2;
        %     plot(rr, ACn);
        %     xlabel('r (\mum)')
        %     ylabel('C(r)');
        %     set(gca,'box','off');
        
    case 2 % --- AUTOCORR FOR CELL OF INTEREST ---
        
        % There is only one point of interest. Get distance between point
        % of interest and all other points
        D = sqrt( (x-x(CellIdx)).^2 + (y-y(CellIdx)).^2 );
        % Discretize
        I = discretize(D, edges);
        
        % --- Perform correlation ---
        
        % For each bin, find all cells in that bin and perform the correlation
        ACp = zeros(1, length(edges)-1)*nan; % preallocate
        ACn = ACp;
        for n=1:length(edges)-1
            I_n = I==n;
            
            theta1 = theta(CellIdx);
            theta2 = theta(I_n);
            
            % Polar autocorrelation. This is mathematically equivalent
            % to an autocorrelation of the unit vectors. This can be shown
            % using the trig identity cos(a-b) = cos(a)cos(b) + sin(a)sin(b).
            ACp_n = sum(cos(theta1-theta2));
            % Nematic autocorrelation.
            ACn_n = sum(cos(2*(theta1-theta2)));
            % Normalize by sum of cos(theta1-theta1) for all points, ie, by the
            % number of points
            ACp(n) = ACp_n / numel(theta2);
            ACn(n) = ACn_n / numel(theta2);
            
        end
        
        % Autocorrelation for shift of 0 is equal to 1
        ACp = [1, ACp];
        ACn = [1, ACn];
end


