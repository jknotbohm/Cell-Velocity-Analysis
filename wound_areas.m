function wound_areas
%
% Calculate area change based on two domain images
%
% If running as a batch, uncomment the statement function at the top and
% comment the clear command
%
% Written by Jacob Notbohm, University of Wisconsin-Madison, 2021
%

% clear;
clc;
close all;

%% --- USER INPUTS ---

% Domain names
d1 = 'domain.tif';
d2 = 'domain_last.tif';

% Name to save data
savename = 'woundarea12.txt';

%% --- GET THE TWO DOMAIN AREAS ---

% Get pixel size Experimental Settings fil
fid = fopen('ExperimentalSettings.txt');
txtcell = cell2mat(textscan(fid,'%f %*[^\n]')); % '%*[^\n]' skips the remainder of each line
pix_size = txtcell(1)*1e6; % Pixel size, um
fclose(fid);

% Get area of wound in first domain image
domain1 = imread(d1);
domain1 = domain1/max(domain1(:));
area1 = sum(1-domain1(:))*pix_size^2; % units: um^2

% Get area of wound in second domain image
domain2 = imread(d2);
domain2 = domain2/max(domain2(:));
area2 = sum(1-domain2(:))*pix_size^2; % units: um^2

%% --- SAVE DATA ---

M = [area1, area2, (area1-area2)/area1]; % Units on areas: um^2
writematrix(M,savename);



