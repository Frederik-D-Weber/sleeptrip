function [rgbcolors labels] = st_epoch_colors(epochs,varargin)

global st_default

colorscheme = st_default.colorscheme; % default
if nargin == 2
    if ~strcmp(varargin{1},'default')
        colorscheme = varargin{1};
    end
end
if ~iscell(epochs)
    epochs = {epochs};
end

lables_order_topdown = {...
    '?','-1',...
    'W','0',...
    'MT','8',...
    'A',...
    'R','5',...
    'S1','N1','1',...
    'S2','N2','2',...
    'S3','N3','SWS','3',...
    'S4','N4','4'};

switch colorscheme
    case 'bright'
    lables_colors_topdown = {...
    [0 0 0],[0 0 0],...
    [0 1 1],[0 1 1],...
    [0 0.5 1],[0 0.5 1],...
    [0.5 0.5 0.5],...
    [0 1 0],[0 1 0],...
    [1 0.75 0],[1 0.75 0],[1 0.75 0],...
    [1 0.5 0],[1 0.5 0],[1 0.5 0],...
    [1 0 0],[1 0 0],[1 0 0],[1 0 0],...
    [1 0 0.75],[1 0 0.75],[1 0 0.75]};
    case 'dark'
    lables_colors_topdown = {...
    [0 0 0],[0 0 0],...
    [0.027, 0.792, 0.792],[0.027, 0.792, 0.792],...
    [0 0.5 1],[0 0.5 1],...
    [0.5 0.5 0.5],...
    [0.007, 0.694, 0.007],[0.007, 0.694, 0.007],...
    [0.878, 0.694, 0],[0.878, 0.694, 0],[0.878, 0.694, 0],...
    [0.898, 0.462, 0.023],[0.898, 0.462, 0.023],[0.898, 0.462, 0.023],...
    [0.784, 0.015, 0.015],[0.784, 0.015, 0.015],[0.784, 0.015, 0.015],[0.784, 0.015, 0.015],...
    [0.709, 0.031, 0.541],[0.709, 0.031, 0.541],[0.709, 0.031, 0.541]}; 
    case 'restless'
        lables_order_topdown = {...
    '?','-1',...
    'W','0',...
    'MT','8',...
    'A',...
    'R','5',...
    'S1','N1','1',...
    'S2','N2','2',...
    'S3','N3','SWS','3',...
    'S4','N4','4'};
    lables_colors_topdown = {...
    [0 0 0],[0 0 0],...
    [0.8 0 0],[0.8 0 0],...
    [0.678 0 0.376],[0.678 0 0.376],...
    [0.5 0.5 0.5],...
    [0.894 0.458 0.027],[0.894 0.458 0.027],...
    [0.027 0.807 0.894],[0.027 0.807 0.894],[0.027 0.807 0.894],...
    [0.019, 0.505, 0.941],[0.019, 0.505, 0.941],[0.019, 0.505, 0.941],...
    [0.019, 0.294, 0.941],[0.019, 0.294, 0.941],[0.019, 0.294, 0.941],[0.019, 0.294, 0.941],...
    [0.015, 0.062, 0.462],[0.015, 0.062, 0.462],[0.015, 0.062, 0.462]};
end
[Lia, Locb] = ismember(epochs,lables_order_topdown);
lables_colors_topdown = lables_colors_topdown(Locb);
rgbcolors = cat(1,lables_colors_topdown{:});

labels = unique(epochs);
[is_in_labels,idx_labels_order] = ismember(labels,lables_order_topdown);
[idx_reorder_topdown,idx_reorder] = sort(idx_labels_order,'descend');
labels = labels(idx_reorder);

end



