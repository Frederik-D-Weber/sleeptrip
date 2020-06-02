function [rgbcolors labels] = st_epoch_colors(epochs)

if ~iscell(epochs)
    epochs = {epochs};
end

lables_order_topdown = {'?','-1','W','0','MT','8','A','R','5','S1','N1','1','S2','N2','2','S3','N3','SWS','3','S4','N4','4'};
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
[Lia, Locb] = ismember(epochs,lables_order_topdown);
lables_colors_topdown = lables_colors_topdown(Locb);
rgbcolors = cat(1,lables_colors_topdown{:});

labels = unique(epochs);
[is_in_labels,idx_labels_order] = ismember(labels,lables_order_topdown);
[idx_reorder_topdown,idx_reorder] = sort(idx_labels_order,'descend');
labels = labels(idx_reorder);

end



