% function find_boundary_scratch_single
%
% Find boundary of a cell layer in an image. This version is for a
% classical scratch wound assay, wherein cells occupy the right and left
% (or top and bottom) sides of the image and migrate to fill the unfilled
% space.
%
% This version also identifies the boundaries of only a single image in a
% time lapse.
%
% This follows some of the recommendations by Treloar & Simpson, Plos One
% e67389, 2013
%
% The script first turns the image into a binary by finding locations in
% the image having gradients larger than a threshold value. Those locations
% are then dilated. Various operations follow to identify the boundary of
% the dilated binary image.
%
% There are 2 important input parameters
%   filter_th -  Threshold used to binarize image based on the gradient.
%   dilate_rad - Radius of circle used to dilate each pixel of the
%                binarized image.
% Before running for every image, the script shows the result of the
% binarization and dilation and asks whether to continue. The parameter you
% will adjust the most is filter_th:
%   - If the binary image is all black, decrease the value of filter_th
%   - If the binary image is all white, increase the value of filter_th
% Sometimes in comparing the binary image to the dilated binary image,
% you'll also find it useful to increase or decrease the value of
% dilate_rad.
%
% This script loads the file 'ExperimentalSettings.txt' and reads the 7th
% row. That row should contain 1 or 2 as follows:
%   1 indicates that the monolayer is a edge/strip geometry with cells on
%   the edge of an image (appropriate for scratch wound assay)
%   2 indicates that the monolayer is an island that is fully within the
%   image.
% This information is needed, because the edge/strip geometry requires some
% extra care in handling pixels at the edge of an image.
%
% Written by Jacob Notbohm, 2015-2021, University of Wisconsin-Madison
%

% This script uses the following functions:
%  smooth2a.m
% This script requires a file called 'ExperimentalSettings.txt.' See readme
% for more information.


clear;
close all;
clc;

mkdir('boundary_plots');

%% --- USER INPUTS ---
% Name of image. Multipage tif is possible. Set to empty arry [] if image
% contains *Pos*.tif
imname = [];
% Threshold for Matlab Sobel filter
filter_th = 6.4e0;
% filter_th = 0.7e3; % Typical value of Hamamatsu Orca Flash: 0.8e3 - 1.2e3
% Radius of circular kernal used to dilate bright spots
dilate_rad = 16;
% Directory and filename to save images with boundary
savenameheader = 'boundary_plots/boundary_';
% Image number in multipage tif to run. If cell image file is not a
% multipage tif, enter 1
image_num = 1;
% Set to [] to make figure window visible
invisible = [];

% Name to save domain
savename = 'domain01.tif';

% State whether to do a check for the filter threshold. Set to empty array
% [] to suppress checking the threshold
filter_th_check = 1;

% Maximum size of holes to remove
% (pixels). These holes appear only on the image boundaries.
N_edge = 1e5; % pix (area)


%% --- LOAD IMAGE ---

% Load image
if isempty(imname)
    F = dir('*Pos*.tif');
    imname = F(1).name;
end

im0 = imread(imname, image_num);
[M, N] = size(im0);


%% --- TEST FILTER THRESHOLD AND DILATE RADIUS ON FIRST IMAGE ---

if ~isempty(filter_th_check)
    
    % Do some minimal filtering
    h = fspecial('gaussian',[3, 3], 0.4);
    im1 = filter2(h,im0);
    % Sometimes a median filter helps to clean up salt and pepper noise, but
    % it's slow and usually isn't needed
    % im1 = medfilt2(im1,[3, 3]);
    % Typically don't use the smoothing option below, but sometimes it helps
    % im1 = smooth2a(im1,5);
    
    % Find edges with sobel filter
    im2 = edge(im1,'sobel',filter_th);
    % Dilate edges
    SE = strel('disk',dilate_rad,0);
    im3 = imdilate(im2, SE);
    
    % Show the two images
    hf = figure;
    pos = get(hf,'Position');
    set(hf,'color','w','position',[0.2*pos(1) 0.8*pos(2) 2.4*pos(3) 0.8*pos(4)]);
    subplot(1,3,1)
    imagesc(im1); colormap gray;
    title('Cell image');
    subplot(1,3,2)
    imagesc(im2); colormap gray;
    title('Binarized image')
    subplot(1,3,3)
    imagesc(im3); colormap gray;
    title('Dilated binary image');
    
    % Ask user if settings are good enough to continue
    answer = questdlg('The dilated binary image should give a rough approximation of the cell monolayer. Continue?',...
        'Continue?', ... % this is the title of the question bar
        'Yes','No','Yes');
    if strcmp(answer,'No')==1
        disp('You selected No. Script canceled.')
        return
    else
        close(hf);
    end
    
end


%% --- FIND BOUNDARIES FOR SPECIFIED TIME POINT ---

% This isn't terribly efficient, because it re-analyzes the image using
% some lines of code above, but for a single image, it's fast

im_k = im0;

% Get rid of some noise in the image. Do only minimal filtering--we
% need as much contrast as possible.
h = fspecial('gaussian',[3, 3], 0.4);
im_kf = filter2(h,im_k);

% --- Optional filtering steps that sometimes help ---
% ! Make sure that the filtering done here matches the filtering done
% in the section titled "TEST FILTER THRESHOLD AND DILATE RADIUS ON
% FIRST IMAGE"

% Sometimes a median filter helps to clean up salt and pepper noise, but
% it's slow and usually isn't needed
% im_k = medfilt2(im_k,[3, 3]);

% Typically don't use the smoothing option below, but sometimes it helps
% im_k = smooth2a(im_k,5);
% ------------------------------------------------------

% Find edges with sobel filter
im_b = edge(im_kf,'sobel',filter_th);

% For testing/debugging: show image
% figure; imagesc(im_b); colormap gray;
%     % For testing/debugging: Find edges with simple threshold
%     im_b = im2bw(uint16(im_k), 200/(2^16-1));
%     im_b = medfilt2(im_b,[5, 5]);
%     figure; imagesc(im_b); colormap gray;

% Dilate edges
SE = strel('disk',dilate_rad,0);
im_b= imdilate(im_b,SE);
%     figure; imagesc(im_1); colormap gray;

% Take inverse to identify "wound" area in middle of image
im_b = ~im_b;

% Median filter. This is slow and generally not needed.
% im_b = medfilt2(im_b,[2*dilate_rad, 2*dilate_rad]);

% --- Get boundaries ---

% If im_b is all 0s, then the domain covers the entire image
if all(im_b(:)==0)
    domain = true(size(im_b));
    B = bwboundaries(domain);
    B_rowcol = B{1};
    xedge = B_rowcol(:,2);
    yedge = B_rowcol(:,1);
else
    
    % Get region with biggest area
    STATS = regionprops(im_b,'Area','PixelIdxList');
    [~, max_idx] = max([STATS.Area]); % Of all connected regions, get max area
    % [~, I] = sort([STATS.Area],'descend'); max_idx = I(1); % Get 2nd largest area
    
    % Get pixels corresponding to region with max area
    idx = STATS(max_idx).PixelIdxList; % linear indices
    B = zeros(size(im_b));
    B(idx)=1;
    B = bwboundaries(B);
    % Take out largest cell in B. This is the set of boundary coordinates
    % for object with largest perimeter.
    [~, max_idx] = max(cellfun('size', B, 1));
    B_rowcol = B{max_idx};
    % Get boundary x and y coordinates
    xedge = B_rowcol(:,2);
    yedge = B_rowcol(:,1);
    
    % Smooth boundary data
    xedge=smooth(xedge);
    yedge=smooth(yedge);
    
    % Set the edge values of 1 to 0. Same for edge
    % values equal to M or N. This is needed so that the poly2mask function
    % includes boundary values. This may have
    % an error of up to 1 pixel at intersections of the image edge and the
    % cell domain edge, but generally we ignore edge data in the stresses,
    % so this shouldn't be a big problem.
    xedge(xedge==1) = 0;
    yedge(yedge==1) = 0;
    xedge(xedge==N) = N+1;
    yedge(yedge==M) = M+1;
    
    % Get all pixels inside the boundary (domain)
    domain = poly2mask(xedge,yedge,M,N);
    
    % Get rid of small junk on the
    % boundaries (smaller than N_edge pixels).
    domain = 1-bwareaopen(1-domain, N_edge);
    
    % Erode region containing cells. Right now, the domain contains area
    % without cells, so dilate it
    domain = imdilate(domain,SE);
    
    % Repeat boundary finding procedure on eroded domain
    B = bwboundaries(domain);
    [~, max_idx] = max(cellfun('size', B, 1));
    B_rowcol = flipud(B{max_idx});
    xedge = B_rowcol(:,2);
    yedge = B_rowcol(:,1);
    xedge = smooth(xedge);
    yedge = smooth(yedge);
    
    % This is to include locations at the edges of an image
    xedge(xedge==1) = 0;
    yedge(yedge==1) = 0;
    xedge(xedge==N) = N+1;
    yedge(yedge==M) = M+1;
    
    domain = poly2mask(xedge,yedge,M,N);
    
    % Take inverse to get area containing cells
    domain = ~domain;
    
end

% Add to cell arrays
xedge_save = xedge;
yedge_save = yedge;

% Make plot
hf = figure;
pos = get(hf,'Position');
set(hf,'color','w','position',[0.2*pos(1) pos(2) pos(3) pos(4)]);
if invisible
    set(hf,'visible','off');
end

imagesc(im_k);
colormap gray; axis equal; axis tight; axis xy;
hold on
plot(xedge,yedge,'g');


%% --- SAVE DOMAIN ---

% Convert to 8-bit integer
domain = uint8(domain)*255;

% Save image
imwrite(domain,savename,'tif'); % saving will overwrite an image with the same name


