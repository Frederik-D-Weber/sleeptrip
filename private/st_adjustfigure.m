function cfg = st_adjustfigure(cfg, fh)

set(fh, 'Name', cfg.title);
title(cfg.title,'Interpreter','none');
figure_width = cfg.figureoutputwidth;     % Width in inches
figure_height = cfg.figureoutputheight;    % Height in inches
pos = get(fh, 'Position');

%set(hhyp, 'Position', [pos(1) pos(2) figure_width*str2num(cfg.figureoutputresolution), figure_height*str2num(cfg.figureoutputresolution)]); %<- Set size
set(fh, 'Position', [pos(1) pos(2) figure_width*100, figure_height*100]); %<- Set size
% Here we preserve the size of the image when we save it.
set(fh,'InvertHardcopy','on');
set(fh,'PaperUnits', cfg.figureoutputunit);

%set(hhyp,'PaperPositionMode','Auto')
set(fh,'PaperSize',[figure_width, figure_height])

papersize = get(fh, 'PaperSize');
left = (papersize(1)- figure_width)/2;
bottom = (papersize(2)- figure_height)/2;
myfiguresize = [left, bottom, figure_width, figure_height];
set(fh,'PaperPosition', myfiguresize);
set(fh,'PaperOrientation', 'portrait');

end