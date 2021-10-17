function h = make_fig(varargin)
% MAKE_FIG.  MAKE FIGURE WINDOW USING NOTBOHM DEFAULTS
% 
% h = make_fig(pos)
% 
% OPTIONAL INPUT:
%       pos is the relative position of the figure window, specified as
%       a ratio of the default figure size. pos is a 4 element vector
%       corresponding to left edge, bottom edge, width, height. If pos is
%       not specified, default values of [0.5 1 0.6 0.6] are used.
% 
% OUTPUT
%       h is a handle to the figure window
% 
% Written by Jacob Notbohm, University of Wisconsin-Madison, 2015
% 

if ~isempty(varargin)
    pos = varargin{1};
else
    pos = [0.5 1 0.6 0.6];
end
    
% Create figure window
h = figure;
% Change background color to white
set(h,'color','w');
% Set position
p = get(h,'position');
set(h,'position',[p(1)*pos(1) p(2)*pos(2) p(3)*pos(3) p(4)*pos(4)]);
% Set fontsize
set(h,'defaultaxesfontsize',11);
set(h,'defaulttextfontsize',11);
% Set font name
set(h,'defaultaxesfontname','arial');
set(h,'defaulttextfontname','arial');
% Set paper position mode for printing
set(h,'paperpositionmode','auto');
% Turn off invert hard copy so that black backgrounds print as black
set(h,'inverthardcopy','off');