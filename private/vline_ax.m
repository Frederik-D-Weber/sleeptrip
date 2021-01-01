function vline_ax(ax, x, varargin)

% VLINE plot a vertical line in the current graph

abc = axis(ax);
x = [x x];
y = abc([3 4]);
if length(varargin)==1
  varargin = {'color', varargin{1}};
end
h = line(x, y, 'Parent', ax);
set(h, varargin{:});
