function [cfg] = st_sleepreport_nin_dutch(cfg, scoring)

% ST_SLEEPREPORT_NIN_DUTCH creates a visual report of sleep data and
% scoring infos in html
%
% Use as
%   [cfg] = st_sleepreport_nin_dutch(cfg, scoring)
%
%
% Optional configuration parameters are:
%  
%  cfg.identifier = a string to name and identify the reported
%                   scoring/recording and give a file name to.
%
% See also ST_SLEEPREPORT, ST_READ_SCORING

% Copyright (C) 2019-, Frederik D. Weber
%
% This file is part of SleepTrip, see http://www.sleeptrip.org
% for the documentation and details.
%
%    SleepTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    SleepTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    SleepTrip is a branch of FieldTrip, see http://www.fieldtriptoolbox.org
%    and adds funtionality to analyse sleep and polysomnographic data.
%    SleepTrip is under the same license conditions as FieldTrip.
%
%    You should have received a copy of the GNU General Public License
%    along with SleepTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id$

ttic = tic;
mtic = memtic;
functionname = getfunctionname();
fprintf([functionname ' function started\n']);

% these are used by the ft_preamble/ft_postamble function and scripts
ft_revision = '$Id$';
ft_nargin   = nargin;
ft_nargout  = nargout;

% do the general setup of the function
st_defaults
ft_preamble init
ft_preamble debug
ft_preamble loadvar data
ft_preamble provenance data
ft_preamble trackconfig

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  return
end

% set defaults

default_identifier = 'report';
if isfield(scoring,'ori')
    if isfield(scoring.ori,'scoringfile')
        [tmppathstr,tmpname,tmpext] = fileparts(scoring.ori.scoringfile);
        if ~isempty(tmpname)
            default_identifier = tmpname;
        end
    end
end


cfg.identifier  = ft_getopt(cfg, 'identifier', default_identifier);
                                    
fprintf([functionname ' function initialized\n']);


fh = figure;


cfg_sd = [];
%cfg_sd.cycle = 'all';
cfg_sd.sleeponsetdef = 'AASM';
res_sd = st_scoringdescriptives(cfg_sd,scoring);

%res_sd.table

cfg_sd = [];
cfg_sd.cycle = 'all';
res_ss = st_sleepcycles(cfg_sd,scoring);

colorscheme = 'restless';

cfg_hp = [];

cfg_hp.colorblocksconnect = 'no';
cfg_hp.plottype = 'colorblocks';
cfg_hp.plotunknown = 'no';
cfg_hp.eventcolorsbystagecolor = 'no';
cfg_hp.eventmaskcolor = [1 1 1];
cfg_hp.colorscheme = colorscheme;


if isfield(scoring,'arousals')
        cfg_hp.eventtimes = {scoring.arousals.start'};
        cfg_hp.eventdurations = {scoring.arousals.duration'};
        %cfg_hp.eventlabels = {'Arousals'};
        cfg_hp.eventlabels = {'Wakker'};
        [wakecolor] = st_epoch_colors('W', colorscheme);
        cfg_hp.eventcolors = wakecolor;
        cfg_hp.eventmask = cfg_hp.eventtimes;
        cfg_hp.eventheight = [1];
        cfg_hp.eventalign = {'center'};
        cfg_hp.eventvalueranges = {[0 1]};
        cfg_hp.eventvalues = {ones(1,numel(cfg_hp.eventtimes{end}))};
        cfg_hp.eventvalueranges_plot = [0];
        cfg_hp.ploteventboundaryticks = [0];
        cfg_hp.offset_event_y = -2.25;
end


if false
        cfg_hp.eventtimes = cat(1,cfg_hp.eventtimes,{(res_ss.table.startepoch(res_ss.table.complete)*scoring.epochlength)'});
        cfg_hp.eventdurations =  cat(1,cfg_hp.eventdurations,{(res_ss.table.durationepochs(res_ss.table.complete)*scoring.epochlength)'});
        cfg_hp.eventlabels = cat(1,cfg_hp.eventlabels,{'Cycles'});
        cfg_hp.eventcolors = cat(1,cfg_hp.eventcolors,[0.5 0.5 0.5]);
        cfg_hp.eventmask = cat(1,cfg_hp.eventmask,zeros(1,numel(cfg_hp.eventtimes{end})));
        cfg_hp.eventheight = cat(1,cfg_hp.eventheight,[1]);
        cfg_hp.eventalign = cat(1,cfg_hp.eventalign, {'bottom'});
        cfg_hp.eventvalueranges = cat(1,cfg_hp.eventvalueranges,{[0 sum(res_ss.table.complete)]});
        cfg_hp.eventvalues = cat(1,cfg_hp.eventvalues,{res_ss.table.cycle(res_ss.table.complete)'});
end

% cfg_hp.eventlabels = {'a','b'};
% cfg_hp.eventtimes = {[30*60 31*60 32*60];[60*60 61*60 62*60]};
% cfg_hp.eventdurations = {[30 30 120];[30 30 120]};
% cfg_hp.eventlabels
% cfg_hp.eventcolors

cfg_hp.plotunknown = 'no';
cfg_hp.plotexcluded = 'no';
cfg_hp.plotsleeponset = 'no';
cfg_hp.plotsleepoffset = 'no';
cfg_hp.plotsleepopon = 'no';
cfg_hp.plotsleepopoff = 'no';
cfg_hp.plotlightsoff = 'no';
cfg_hp.plotlightson = 'no';

cfg_hp.timeunitdisplay  = 'hours';
cfg_hp.timeticksdiff = 60;
cfg_hp.plotlegend = 'no';
cfg_hp.timeunitdisplay = 'time';

cfg_hp.relabelstages = {'Wakker', 'REM', 'Lichte slaap (N1)', 'Diepe slaap (N2)', 'Zeer diepe slaap (N3)'};
%cfg_hp.plottype = 'colorbar';
%cfg_hp.relabelstages = {'Wakker'};

cfg_hp.figureoutputfile       = [cfg.identifier '.png'];
cfg_hp.figureoutputformat     = 'png';

cfg_hp.timestamp              = 'no';
cfg_hp.folderstructure        = 'no';
cfg_hp.xlabel                 = 'Tijd';
cfg_hp.ylabel                 = 'Slaap stadium';
cfg_hp.title                  = datestr(scoring.startdatetime,'dd-mmm-yyyy HH:MM:SS');

[fh_hp axh_hp] = st_hypnoplot(cfg_hp, scoring);

close(fh_hp)

hypnogram_file_path = cfg_hp.figureoutputfile;




hr_min_str = @(mins) [sprintf('%d',fix(mins/60)) ' uur en ' sprintf('%d',fix(rem(mins,60))) ' minuten'];   
percent_str = @(perc) [sprintf('%5.1f',round(perc,1)) '%'];   
density_per_hour_str = @(dens_hour) [sprintf('%5.1f',round(dens_hour,1)) ' per uur'];   

%'  border:1px solid black;',char(10),...

report_html_text = ...
['',...
'<!DOCTYPE html>',char(10),...
'<html>',char(10),...
'<style>',char(10),...
'table, th, td {',char(10),...
'  text-align: left;',char(10),...
'}',char(10),...
'</style>',char(10),...
'<body id="body">',char(10),...
'<table style="width:100%">',char(10),...
'  <tr>',char(10),...
'    <th>Uw nacht in cijfers</th>',char(10),...
'    <th>In tijd</th>',char(10),...
'    <th>In percentage</th>',char(10),...
'  </tr>',char(10),...
'  <tr>',char(10),...
'    <td>Wakker voor in slaap vallen</td>',char(10),...
'    <td> ' hr_min_str(res_sd.table.sleep_onset_delay_min) '</td>',char(10),...
'    <td></td>',char(10),...
'  </tr>',char(10),...
'  <tr>',char(10),...
'    <td>Slaap periode inclusief tussendoor wakker liggen</td>',char(10),...
'    <td>' hr_min_str(res_sd.table.total_sleep_period_duration_min) '</td>',char(10),...
'    <td></td>',char(10),...
'  </tr>',char(10),...
'    <tr>',char(10),...
'    <td>Periodes tussendoor wakker liggen</td>',char(10),...
'    <td>samen ' hr_min_str(res_sd.table.Wake_after_sleep_onset_of_sleep_period_min) '</td>',char(10),...
'    <td>' percent_str(res_sd.table.Wake_after_sleep_onset_perc_of_sleep_period) '</td>',char(10),...
'  </tr>',char(10),...
'  </tr>',char(10),...
'    <tr>',char(10),...
'    <td>Lichte slaap (N1)</td>',char(10),...
'    <td>' hr_min_str(res_sd.table.N1_of_sleep_period_min) '</td>',char(10),...
'    <td>' percent_str(res_sd.table.N1_perc_of_sleep_period) '</td>',char(10),...
'  </tr>',char(10),...
'  </tr>',char(10),...
'    <tr>',char(10),...
'    <td>Diepe slaap (N2)</td>',char(10),...
'    <td>' hr_min_str(res_sd.table.N2_of_sleep_period_min) '</td>',char(10),...
'    <td>' percent_str(res_sd.table.N2_perc_of_sleep_period) '</td>',char(10),...
'  </tr>',char(10),...
'  </tr>',char(10),...
'    <tr>',char(10),...
'    <td>Zeer diepe slaap (N3)</td>',char(10),...
'    <td>' hr_min_str(res_sd.table.N3_of_sleep_period_min) '</td>',char(10),...
'    <td>' percent_str(res_sd.table.N3_perc_of_sleep_period) '</td>',char(10),...
'  </tr>',char(10),...
'  </tr>',char(10),...
'    <tr>',char(10),...
'    <td>REM slaap</td>',char(10),...
'    <td>' hr_min_str(res_sd.table.R_of_sleep_period_min) '</td>',char(10),...
'    <td>' percent_str(res_sd.table.R_perc_of_sleep_period) '</td>',char(10),...
'  </tr>',char(10),...
'  </tr>',char(10),...
'    <tr>',char(10),...
'    <td>Onderbrekingen uit REM slaap</td>',char(10),...
'    <td>' sprintf('%d',res_sd.table.interrupts_W_MT_NR_number_in_all_R_episodes)  ' keer</td>',char(10),...
'    <td>' percent_str((res_sd.table.interrupts_W_MT_NR_overall_durations_min_in_all_R_episodes)/ res_sd.table.R_of_sleep_period_min) ' van de REM</td>',char(10),...
'  </tr>',char(10),...
'  </tr>',char(10),...
'    <tr>',char(10),...
'    <td>REM onrust</td>',char(10),...
'    <td>' density_per_hour_str(res_sd.table.REM_fragmentation_density_per_hour_in_R_epsds) ' van de REM</td>',char(10),...
'    <td></td>',char(10),...
'  </tr>',char(10),...
'</table>',char(10),...
'',...
'<figure>',char(10),...
'  <img src="' hypnogram_file_path '" alt="Hypnogram" style="width:100%">',char(10),...
'</figure>',char(10),...
'</body>',char(10),...
'</html>',char(10),...
''];

% '<script src="https://cdnjs.cloudflare.com/ajax/libs/jspdf/1.3.2/jspdf.min.js"></script>',char(10),...
% '<script>',char(10),...
% '    function downloadPDF() {',char(10),...
% '      var doc = new jsPDF();',char(10),...
% '',char(10),...
% '      doc.fromHTML(document.getElementById(''body''), 15, 15, {',char(10),...
% '        ''width'': 170, ',char(10),...
% '      });',char(10),...
% '      doc.save();',char(10),...
% '    }',char(10),...
% '</script>',char(10),...
% '<button onclick="downloadPDF()">PDF</button>',char(10),...

%'<script src="https://unpkg.com/jspdf@latest/dist/jspdf.umd.min.js"></script>',char(10),...
%'<button id="button">PDF</button>',char(10),...
%'<script src="html2pdf.bundle.min.js"></script>',char(10),...

% '<a href="" download="' [cfg.identifier '.pdf'] '">Download as PDF</a>',char(10),...
% '  <figcaption>' 'slaap datum en tijd' datestr(scoring.startdatetime) '</figcaption>',char(10),...

% '<script>',char(10),...
% 'const btn = document.getElementById("button");',char(10),...
% '',char(10),...
% 'btn.addEventListener("click", function(){',char(10),...
% 'var element = document.getElementById(''body'');',char(10),...
% 'html2pdf().from(element).save(''' [cfg.identifier '.pdf'] '.pdf'');',char(10),...
%'});',char(10),...

% https://github.com/eKoopmans/html2pdf.js

file_html = fopen([cfg.identifier '.html'],'w');
fprintf(file_html,'%s',report_html_text);
fclose(file_html);

% 
% subplot_rows = 2;
% subplot_columns = 1;
% subplot(subplot_rows,subplot_columns,1);
% 
% 
% 
% %figure
% sd_names = {'SOL' 'SOL' 'WASO' 'N1' 'N2' 'N3' 'REM'};
% 
% sd_values = [-res_sd.table.sleep_onset_delay_min res_sd.table.sleep_onset_delay_min res_sd.table.Wake_after_sleep_onset_of_sleep_period_min res_sd.table.N1_of_sleep_period_min res_sd.table.N2_of_sleep_period_min res_sd.table.N3_of_sleep_period_min res_sd.table.R_of_sleep_period_min];
% sd_values = sd_values/60;
% 
% bh = barh([sd_values;sd_values],'stacked');
% axh_b = gca;
% set(axh_b,'YLim',[0.5 1.5])
% hold(axh_b,'on');
% 
% bh(1).FaceColor = [0.8 0.8 0.8];
% bh(2).FaceColor = [0.8 0.8 0.8];
% bh(3).FaceColor = st_epoch_colors('W', colorscheme);
% bh(4).FaceColor = st_epoch_colors('N1', colorscheme);
% bh(5).FaceColor = st_epoch_colors('N2', colorscheme);
% bh(6).FaceColor = st_epoch_colors('N3', colorscheme);
% bh(7).FaceColor = st_epoch_colors('R', colorscheme);
% 
% tempoffset_x = cumsum([0 sd_values]);
% for iSD = 2:numel(sd_names)
%     text(tempoffset_x(iSD),1.5-(iSD-1)*(1/(numel(sd_names))),...
%         sprintf([sd_names{iSD} ': %-5.2f h %-5.2f%%'],sd_values(iSD),100*sd_values(iSD)/(res_sd.table.total_sleep_period_duration_min/60)),...
%         'HorizontalAlignment','left',...
%         'FontSize',12)
% end
% 
% hold(axh_b,'off');
% 
% 
% subplot(subplot_rows,subplot_columns,2);
% 
% figure
% 
% axh_t = gca;
% text(1,1.1,sprintf('Slaap efficiëntie: %-5.2f%%',res_sd.table.sleep_eff_total_sleep_dur_in_sleep_prd_perc_of_sleep_prd),...
%         'HorizontalAlignment','left',...
%         'FontSize',12)
% 
% if isfield(scoring,'arousals')
%     text(1,1,sprintf('Ontwaken: %d per nacht gedurende in totaal %-5.2f minuten',res_sd.table.arousals_number, sum(scoring.arousals.duration)/60),...
%         'HorizontalAlignment','left',...
%         'FontSize',12)
% end
% 
% 
% set(gca,'YLim',[0.9 1.3])
% set(gca,'XLim',[0.5 4])
% 
% 

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble trackconfig
% ft_postamble previous data
% data = datanew
% ft_postamble provenance data
% ft_postamble history    data
% ft_postamble savevar    data


fprintf([functionname ' function finished\n']);
toc(ttic)
memtoc(mtic)
end
