function [AC, G] = autocorr_nogrid(varargin)
%AUTOCORR_NOGRID Autocorrelation of data not on a uniform grid
%
% AC = autocorr_nogrid(x, y, Tx, Ty, Npts, spacing, CellIdx)
% or
% [AC, G] = autocorr_nogrid(x, y, Tx, Ty, Npts, spacing, CellIdx)
%
% This function performs the correlation by binning data and computing the
% autocorrelation of all data within each bin.
%
% This function runs autocorrelation on a vector quantity. Note: Don't use
% this function on unit vectors (which will give an autocorrelation of
% angles). Instead, use autocorr_angle_nogrid.m.
%
% INPUTS
% x, y      vectors identifying points
% Tx, Ty    vectors identifying variable to autocorrelate
% Npts      Number of points at which to compute the autocorrelation. This
%           gives the length of the autocorrelation vector
% spacing   Spacing between the bins in the autocorrelation vector
% CellIdx   (Optional parameter) Scalar corresponding to index of the cell
%           for which to compute the autocorrelation. If no value is input,
%           the autocorrelation is computed and averaged for all cells,
%           which is likely the more common usage of this function.
%
% OUTPUT
% AC: Normalized autocorrelation
%
% OPTIONAL OUTPUT
% G:  Number of data points--used to get pair correlation function
%
% Note: The output AC should start with a value of 1, ie, AC(1)==1.
% Subsequent values of AC should not be NaN; if they are NaN, then no
% points exist within the requested spacing--increase the spacing.
%
% Written by Jacob Notbohm, University of Wisconsin-Madison 2020-2021
%

% Check outputs
if nargout>2
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

% For autocorrelation, need to subtract mean values.
% ! If you input unit vectors, you don't want to subtract off mean values.
% ! Hence you should not use this function if unit vectors are input.
Tx2 = Tx - mean(Tx(:));
Ty2 = Ty - mean(Ty(:));

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
        
        % Since D is symmetric, in principle, we only need to consider the
        % upper triangular portion. This could avoid double counting cells and will
        % cut computation time in half. However, using the upper triangular portion
        % will affet how the normalization of the autocorrelation is performed--if
        % the upper triangular portion is kept, the normalization based on I1 (see
        % I1 defined below) will be missing some data points. Keeping the entire
        % matrix D will avoid this issue, and I1 will contain all data, enabling us
        % to use I1 in the normalization.
        
        % --- Perform correlation ---
        
        I = discretize(D, edges);
        % nan values occur for cells that are a distance beyond the max of edges
        
        % For each bin, find all cells in that bin and perform the correlation
        AC = zeros(1, length(edges)-1); % preallocate
        G = AC;
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
                AC1 = sum(Tx2(I1).*Tx2(I2) + Ty2(I1).*Ty2(I2));
                
                % Denominator
                % The typical option for the denominator is below. Note,
                % however, that if this option is chosen, you should NOT pick the
                % upper triangular region of D; you need to use the full matrix.
                AC2 = sum(Tx2(I1).^2 + Ty2(I1).^2);
                
                % normalized autocorrelation
                AC(n) = AC1/AC2;
                
                % For two point correlation function
                G(n) = length(I1);
                
            end
        end
        
        % Autocorrelation for shift of 0. This is a correlation of all points with
        % themselves and then normalized according to equations above.
        AC1 = sum(Tx2.*Tx2 + Ty2.*Ty2);
        AC2 = sqrt(sum(Tx2.^2))*sqrt(sum(Tx2.^2)) + sqrt(sum(Ty2.^2))*sqrt(sum(Ty2.^2));
        
        AC = [AC1/AC2, AC];
        
        
        
        %     % Plot for checking/debugging
        %     figure;
        %     rr = (0:1:Npts-1)*spacing + spacing/2;
        %     plot(rr, AC);
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
        AC = zeros(1, length(edges)-1); % preallocate
        G = AC;
        for n=1:length(edges)-1
            I_n = find(I==n); % Column vector
            
            % Perform autocorrelation
            % numerator
            AC1 = sum(Tx2(CellIdx)*Tx2(I_n) + Ty2(CellIdx)*Ty2(I_n));
            
            % Denominator isn't cleanly defined...below is a reasonable
            % estimate but nor perfectly accurate, making this version of
            % the script not fully rigorous
            AC2 = sum( sqrt(Tx2(CellIdx)*Tx2(I_n)) + sqrt(Ty2(CellIdx)*Ty2(I_n)) );
            
            % normalized autocorrelation
            AC(n) = AC1/AC2;
            
            % For two point correlation function
            G(n) = length(I_n);
            
        end
        
        % Add autocorrelation for shift of 0
        AC = [1, AC];
end


% For two point correlation function, data for spacing of 0 is nan
if nargout==2
    G = [nan, G];
end
