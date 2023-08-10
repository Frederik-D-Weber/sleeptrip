function fh=st_plot_data_quality(cfg)

% Copyright (C) 2022-, Roy Cox, Frederik D. Weber
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
%ft_preamble init
% ft_preamble debug
% ft_preamble loadvar data
% ft_preamble provenance data
% ft_preamble trackconfig

% the ft_abort variable is set to true or false in ft_preamble_init
% if ft_abort
%   return
% end

%-----core function start---
%---input checks and defaults----
ft_checkconfig(cfg,'required',{'grid','artifact_summary'});
cfg.title  = ft_getopt(cfg, 'title', '');
fig_title=cfg.title;
cfg.style  = ft_getopt(cfg, 'style', 'full');

%get the grid
cfg_grid=cfg.grid;

%check if plotting style is permitted
detector_labels=cfg_grid.label;
if ~ismember(cfg.style,[detector_labels 'full' 'basic'])
    cfg.style='full'; %default to full
end
plot_style=cfg.style;

fprintf([functionname ' function initialized\n'])


hasScoring=false;
if isfield(cfg,'scorings')
    scoring=cfg.scorings.scoring_artifact_level;
    hasScoring=true;
end



numChan=cfg_grid.channel_number;
numSegment=cfg_grid.segment_number;
segmentLength=cfg_grid.segment_length;
dataDuration=(numSegment*segmentLength)/60;
chanLabel=cfg.elec.label;

%----define some colors
white=[1 1 1];
black=[0 0 0];
dark_gray=[0.5 0.5 0.5];
light_gray=[0.7 0.7 0.7];
red=[255 50 50]/255;
dark_red=[200 0 0]/255; %other red
beige=[245 245 220]/255;
brown=[218,165,32]/255;
dark_brown=[139,69,19]/255;
light_brown=[222,184,135]/255;
deep_pink=[255,105,180]/255;
hot_pink=[255,105,180]/255;
magenta=[1 0 1];
medium_purple=[147,112,219]/255;
indigo=[75,0,130]/255;
yellow=[1 1 0];
gold=[255,215,0]/255;
orange=[255,150,150]/255;
dark_green=[0 100 0]/255;

%basic colors
figure_col=white;

%assign colors to artifact types and build colormaps
default_col=white;
artifact_type_labels={'basic','spatial exp.','temporal exp.'};
artifact_cols=[dark_gray;light_brown; brown];
reject_label={'REJECT'};
reject_col=red;


artifact_type_reject_labels=[artifact_type_labels reject_label];
colormap_bar=[artifact_cols; reject_col];
colormap_grid=[default_col; artifact_cols; reject_col];



%topography colormap
startColTopo=white;
endColTopo=black;
colormap_topo=[linspace(startColTopo(1),endColTopo(1))' linspace(startColTopo(2),endColTopo(2))' linspace(startColTopo(3),endColTopo(3))'];


xax_col=black;
yax_col=black;

%----setup figure--
fh=figure('color',figure_col,'Units','centimeters','Position',[2 2 50 25]);

fig_rows=9;
fig_cols=18;
subplot_layout=reshape(1:fig_rows*fig_cols,[fig_cols fig_rows])';


table_rows=1:fig_rows;
hyp_rows=1;
grid_rows=hyp_rows(end)+1:7; %2-7
seg_rows=grid_rows(end)+1:fig_rows;

table_cols=1:3;
grid_cols=4:15;

%summary table
ax_table=subplot(fig_rows,fig_cols,reshape(subplot_layout(table_rows,table_cols),[],1));
%hypnogram
ax_hyp=subplot(fig_rows,fig_cols,reshape(subplot_layout(hyp_rows,grid_cols),[],1));
%hypnogram legend plus sleep architecture
ax_leg=subplot(fig_rows,fig_cols,reshape(subplot_layout(hyp_rows,grid_cols(end)+1:fig_cols),[],1));
%grid
ax_grid=subplot(fig_rows,fig_cols,reshape(subplot_layout(grid_rows,grid_cols),[],1));
%channel summary
ax_chan=subplot(fig_rows,fig_cols,reshape(subplot_layout(grid_rows,grid_cols(end)+1:fig_cols),[],1));
%segment summary
ax_seg=subplot(fig_rows,fig_cols,reshape(subplot_layout(seg_rows,grid_cols),[],1));
%topography
ax_topo=subplot(fig_rows,fig_cols,reshape(subplot_layout(seg_rows,grid_cols(end)+1:fig_cols),[],1));



%plotting details
fsize=10; %general fontsize
fsize_legend=8;
fsize_small=7;
fsize_big=12;

maxChanLabels=30; %max channel labels to plot

reject_alpha=0.5;
ax_width=1;
%%

%precalculate some things

bad_chan_inds=find(mean(cfg.grid.artifact_grid_segment_expanded,2)==1);

reject_grid=double(cfg_grid.reject_grid);
rejectprop=mean(reject_grid(1,:),2);
badchannelthresh_adjusted=cfg.badchannelthresh*(1-rejectprop)+rejectprop;
%%
axes(ax_table);



ycoor_start=0;

ycoor_int=0.4;ycoor_int_small=0.2;ycoor_int_big=0.6;
xleft=0;xright=3;

recordingDetails={'plot style' plot_style;...
    'RECORDING DETAILS' '';...
    'data length' [num2str(dataDuration,'%.1f') ' min'];...
    'segments' num2str(numSegment,'%i');...
    'segment length' [num2str(segmentLength,'%i') ' s'];...
    'channels' num2str(numChan,'%i')};


summaryTable=cfg.artifact_summary;

%find all detector values (in convenient order)
rowNames=summaryTable.Properties.RowNames;
[~, detectInds]=ismember(['REJECT','REPAIR','BASIC',detector_labels,'spatial expansion','temporal expansion','ANY'],rowNames);
detectNames=rowNames(detectInds);
detectVals=num2cell(summaryTable{detectInds,'artifact_data_global_perc'});

%rename some
detectNames(strcmp(detectNames,'BASIC'))={'basic'};
detectNames(strcmp(detectNames,'spatial expansion'))={'spatial exp.'};
detectNames(strcmp(detectNames,'temporal expansion'))={'temporal exp.'};

detectStr= cellfun(@(X) [num2str(X,'%.1f') '%'],detectVals,'UniformOutput',false);

artifactDetails=[detectNames detectStr];

%set up all string info
allDetails=[recordingDetails; artifactDetails];
%add before REJECT:
insert_ind=find(strcmp(allDetails(:,1),'REJECT'));
allDetails=[allDetails(1:insert_ind-1,:) ; {'ARTIFACT SUMMARY' ''} ;allDetails(insert_ind:end,:)];
%add before REPAIR:
insert_ind=find(strcmp(allDetails(:,1),'REPAIR'));
allDetails=[allDetails(1:insert_ind-1,:) ; {'threshold' [num2str(100*cfg.segmentrejectthresh,'%.1f') '%']} ;allDetails(insert_ind:end,:)];
%add before basic:
insert_ind=find(strcmp(allDetails(:,1),'basic'));
allDetails=[allDetails(1:insert_ind-1,:) ; {'bad channels' [num2str(length(bad_chan_inds),'%i') ' (' num2str(100*length(bad_chan_inds)/numChan,'%.1f') '%)' ]} ;allDetails(insert_ind:end,:)];
%add before basic:
insert_ind=find(strcmp(allDetails(:,1),'basic'));
allDetails=[allDetails(1:insert_ind-1,:) ; {'ARTIFACT TYPES' ''} ;allDetails(insert_ind:end,:)];
%add before temporal exp:
insert_ind=find(strcmp(allDetails(:,1),'temporal exp.'));
allDetails=[allDetails(1:insert_ind-1,:) ; {'threshold' [num2str(100*cfg.channelexpandthresh,'%.1f') '%']} ;allDetails(insert_ind:end,:)];
%add before ANY:
insert_ind=find(strcmp(allDetails(:,1),'ANY'));
allDetails=[allDetails(1:insert_ind-1,:) ; {'threshold' [num2str(100*cfg.badchannelthresh,'%.1f') '%']} ;allDetails(insert_ind:end,:)];
% %add before ANY:
% insert_ind=find(strcmp(allDetails(:,1),'ANY'));
% allDetails=[allDetails(1:insert_ind-1,:) ; {'threshold (true)' [num2str(100*badchannelthresh_adjusted,'%.1f') '%']} ;allDetails(insert_ind:end,:)];



stringNames=allDetails(:,1);

%fontsize: standard unless individual detector
detectFontSize=num2cell(ones(size(allDetails,1),1)*fsize);
detectFontSize(ismember(stringNames,[detector_labels 'threshold' 'threshold (raw)' 'threshold (true)' 'bad channels']))={fsize_small};
detectFontSize(ismember(stringNames,{'REJECT','REPAIR'}))={fsize};


%color: black unless reject or repair
detectFontColor=repmat({'k'},size(allDetails,1),1);
detectFontColor(ismember(stringNames,{'REJECT'}))={red};
detectFontColor(ismember(stringNames,{'REPAIR'}))={magenta};
detectFontColor(ismember(stringNames,{'basic'}))={artifact_cols(1,:)};
detectFontColor(ismember(stringNames,{'spatial exp.'}))={artifact_cols(2,:)};
detectFontColor(ismember(stringNames,{'temporal exp.'}))={artifact_cols(3,:)};


%%line interval (before item): standard unless individual detector
detectLineInt=num2cell(ones(size(allDetails,1),1)*ycoor_int); %standard
detectLineInt(ismember(stringNames,[detector_labels 'threshold' 'threshold (raw)' 'threshold (true)' 'bad channels']))={ycoor_int_small}; %smaller for individual detectors
detectLineInt(ismember(stringNames,{'RECORDING DETAILS','ARTIFACT SUMMARY','ARTIFACT TYPES'}))={ycoor_int_big}; %larger for headers
detectLineInt(1)={1};%except first item

%bold/normal
detectFontweight=repmat({'normal'},size(allDetails,1),1);
detectFontweight(ismember(stringNames,{'RECORDING DETAILS','ARTIFACT SUMMARY','ARTIFACT TYPES','basic','spatial exp.' ,'temporal exp.' ,'ANY','REJECT','REPAIR' }))={'bold'};

%combine all properties and loop to print
allTextProps=[allDetails detectFontSize detectFontColor detectLineInt detectFontweight];

for text_i=1:size(allTextProps,1)
    ycoor_start=ycoor_start+allTextProps{text_i,5};
    text(xleft, ycoor_start,allTextProps{text_i,1},'FontSize',allTextProps{text_i,3},'HorizontalAlignment','left','Color',allTextProps{text_i,4},'Fontweight',allTextProps{text_i,6})
    text(xright, ycoor_start,allTextProps{text_i,2},'FontSize',allTextProps{text_i,3},'HorizontalAlignment','right')

end


g=gca;
set(g,'TickDir','none',...
    'XTick',[],...
    'YTick',[],...
    'YDir','reverse',...
    'box','off','Visible','off')

xlim([0 5])
ylim([0 10])




%%
%----HYPNOGRAM---

if ~hasScoring %no scoring

    %plot "no scoring" in place of hypnogram
    axes(ax_hyp);

    text(0,0,'no scoring','FontSize',fsize_big,'HorizontalAlignment','center','Color','k')

    g=gca;
    set(g,'TickDir','none',...
        'XTick',[],...
        'YTick',[],...
        'box','off','Visible','off')

    xlim([-1 1]);ylim([-1 1]);

    %plot nothing for legend
    axes(ax_leg);

    g=gca;
    set(g,'TickDir','none',...
        'XTick',[],...
        'YTick',[],...
        'box','off','Visible','off')
elseif hasScoring

    axes(ax_hyp);

    %get simple sleep architecture
    tmp_cfg=[];
    tmp_cfg.sleeponsetdef='XR'; %any stage
    descriptives=st_scoringdescriptives(tmp_cfg,scoring);

    stage_labels = unique(scoring.epochs);
    stage_label_colors = st_epoch_colors(stage_labels, 'restless');

    %set up a numeric vector for plotting
    scoring_image=zeros(size(scoring.epochs));
    stage_minutes=nan(size(stage_labels));
    for label_i=1:length(stage_labels)
        curr_stage_label=stage_labels{label_i};

        scoring_image(strcmp(scoring.epochs,curr_stage_label))=label_i;

        if ismember(curr_stage_label,{'N1','N2','N3','R'})
            var_name_mins=[curr_stage_label '_all_scoring_min'];
        elseif strcmp(curr_stage_label,'W')
            var_name_mins='Wake_all_scoring_min';
        elseif strcmp(curr_stage_label,'?')
            var_name_mins='Unknown_all_scoring_min';
        else
            continue

        end
        stage_minutes(label_i)=descriptives.table.(var_name_mins);

    end
    scoring_duration=descriptives.table.scoring_duration_min;

    stage_percentage=stage_minutes/scoring_duration*100;

    imagesc(scoring_image)
    g=gca;


    colormap(g,stage_label_colors);

    set(g,'TickDir','none',...
        'XTick',[],...
        'YTick',[],...
        'box','off','Visible','off')

    %-------HYPNOGRAM "LEGEND" WITH SLEEP ARCHITECTURE----
    axes(ax_leg);


    %make the labels
    legend_labels=[stage_labels;...
        cellfun(@(X) [num2str(X,'%.1f') ' min'],num2cell(stage_minutes),'UniformOutput',false);...
        cellfun(@(X) [num2str(X,'%.1f') '%'],num2cell(stage_percentage),'UniformOutput',false)];

    legend_labels=mat2cell(legend_labels,3,ones(1,length(stage_labels)));


    text_y_coor=1;

    startPoints=1:length(stage_labels);
    xDat=[startPoints' startPoints'+0.9];
    yDat=ones(size(xDat));

    plot(xDat',yDat','LineWidth',40)

    xlim([min(xDat(:)) max(xDat(:))]);


    g=gca;

    set(g,'TickDir','none',...
        'XTick',[],...
        'YTick',[],...
        'box','off',...
        'ColorOrder',stage_label_colors,...
        'Visible','off')

    text(mean(xDat,2),repmat(text_y_coor,[1 length(stage_labels)]),legend_labels,...
        'FontSize',fsize_legend,'HorizontalAlignment','center','FontWeight','bold','Color','w');

    ylim([0.5 1.5])

end


%%
%-------------GRID------------
axes(ax_grid);
hold on

%build artifact_grid_plot (non-verlapping)
if ismember(plot_style,{'basic' 'full'})

    %"basic" artifact type (union of all individual ones)
    artifact_grid_plot=double(cfg_grid.artifact_grid_merged);    %0 (clean) and 1 (basic artifact)

    %store for topography
    artifact_data_topo=artifact_grid_plot;

    if strcmp(plot_style,'full')

        %extend to spatial and temporal expansion
        artifact_grid_plot(cfg_grid.artifact_grid_channel_expansion==1)=2; %2)channel expanded
        artifact_grid_plot(cfg_grid.artifact_grid_segment_expansion==1)=3;%3)segment expanded (entire channel)

        %set up legend
        grid_legend_labels=artifact_type_reject_labels;
        legendOrder=[4 1:3];

        %values in artifact_grid_plot to use for downstream bar plots
        evalGridVals=1:3;

    elseif strcmp(plot_style,'basic')

        %set up legend
        grid_legend_labels=artifact_type_reject_labels(1);
        legendOrder=1;

        %values in artifact_grid_plot to use for downstream bar plots
        evalGridVals=1;

    end
else
    artifact_grid_plot=double(squeeze(cfg_grid.artifact_grid_by_type(strcmp(detector_labels,plot_style),:,:))); %0 (clean) and 1 (basic artifact)

    artifact_data_topo=artifact_grid_plot;

    %set up legend
    grid_legend_labels={plot_style};
    legendOrder=1;

    %values in artifact_grid_plot to use for downstream bar plots
    evalGridVals=1;

end


%plot different artifact types (scale 0-4, regardless of actual max value)
imagesc(1:numSegment, 1:numChan,artifact_grid_plot,[0 4])

%assign colormap
g=gca;
colormap(g,colormap_grid);

%plot transparent rejection grid on top (scale 0-4)
if strcmp(cfg.style,'full')


    %reject_grid
    reject_grid_plot=reject_grid;
    reject_grid_plot(reject_grid_plot==0)=nan; %convert 0 to nan so that 0 corresponds to no color
    reject_grid_plot(reject_grid_plot==1)=4;

    %transparency grid (fully transparent, unless reject==1)
    reject_alpha_grid=zeros(size(reject_grid_plot));
    reject_alpha_grid(reject_grid==1)=reject_alpha;

    %add rejection grid with transparency mask
    imagesc(1:numSegment, 1:numChan,reject_grid_plot,'AlphaData',reject_alpha_grid,[0 4])
end

%customize plot appearance
xticks=[];
label_interval=ceil(numChan/maxChanLabels);
chan_ticks=1:label_interval:numChan;
chan_ticklabels=chanLabel(chan_ticks);


g=gca;
set(g,'TickDir','out',...
    'XTick',xticks,...
    'YTick',chan_ticks,'YTickLabel',chan_ticklabels,...
    'FontSize',fsize,'box','off',...
    'YDir','reverse',...%'TickLength',[0.01 0.01],...
    'XColor',xax_col,'YColor',yax_col)
g.XAxis.LineWidth=ax_width;
g.YAxis.LineWidth=ax_width;
xlim([0.5 numSegment+0.5])
ylim([0.5 numChan+0.5])
ylabel('channel')

%%
%----------------CHANNEL-WISE---------
axes(ax_chan);
hold on

chanHealth=cell2mat(arrayfun(@(X) 100*mean(artifact_grid_plot==X,2),evalGridVals,'UniformOutput',false));
b=barh(1:numChan,flipud(chanHealth),'stacked','FaceColor','flat','EdgeColor','none','BarWidth',1);
for b_i=1:length(b)
    b(b_i).CData = colormap_bar(b_i,:);
end

if strcmp(cfg.style,'full')
    b(end+1)=barh(1:numChan,100*mean(reject_grid,2),'FaceColor',colormap_bar(end,:),'FaceAlpha',reject_alpha,'EdgeColor','none','BarWidth',1);

    %add threshold used for bad channel expansion (adjusted)

%     H=vline(100*badchannelthresh_adjusted);
%     H.LineStyle='--';
%     H.LineWidth=1;
%     H.Color=dark_green;
end

xlim([0 100])
ylim([0.5 numChan+0.5])
xlabel('segments (%)')
g=gca;
set(g,'TickDir','out',...
    'XAxisLocation','bottom',...
    'XTick',0:25:100,'XGrid','on',...
    'YTick',[],...%'TickLength',[0.01 0.01],...
    'FontSize',fsize,'box','off','XColor',xax_col,'YColor',yax_col,'LineWidth',1)
g.XAxis.LineWidth=ax_width;
g.YAxis.LineWidth=ax_width;


legend(b(legendOrder),grid_legend_labels(legendOrder),'box','on','FontSize',fsize,'Location','northeast','FontWeight','normal')
%%
%-----------SEGMENT-WISE--------
axes(ax_seg);
hold on

segmentHealth=cell2mat(arrayfun(@(X) 100*mean(artifact_grid_plot==X,1),evalGridVals,'UniformOutput',false)')';
b=bar(1:numSegment,segmentHealth,'stacked','FaceColor','flat','EdgeColor','none','BarWidth',1);
for b_i=1:length(b)
    b(b_i).CData = colormap_bar(b_i,:);
end

if strcmp(cfg.style,'full')
    b(end+1)=bar(1:numSegment,100*mean(reject_grid,1),'FaceColor',colormap_bar(end,:),'FaceAlpha',reject_alpha,'EdgeColor','none','BarWidth',1);

%     %add threshold used for rejection
%     V=hline(100*cfg.segmentrejectthresh);
%     V.LineStyle='--';
%     V.LineWidth=1;
%     V.Color=dark_green;
end

%calculate time vector for ticks
time_vect_sec=segmentLength*((1:numSegment)-1);

tick_interval_sec=60*30;
duration_vect=seconds(time_vect_sec);
duration_vect.Format='hh:mm';

num_intervals=max(floor(time_vect_sec/tick_interval_sec));
xticks=dsearchn(time_vect_sec',[0:num_intervals]'*tick_interval_sec);
xticklabels=char(duration_vect(xticks));

xlim([0.5 numSegment+0.5])
ylim([0 100])
ylabel('channels (%)')
xlabel('elapsed time (hh:mm)')
g=gca;
set(g,'TickDir','out',...
    'XTick',xticks,'XTickLabel',xticklabels,...
    'YTick',0:25:100,'YGrid','on',...
    'FontSize',fsize,'box','off',...%'TickLength',[0.01 0.01],...
    'XColor',xax_col,'YColor',yax_col)
g.XAxis.LineWidth=ax_width;
g.YAxis.LineWidth=ax_width;


%%
%---------------TOPOGRAPHY-------------------
axes(ax_topo);

%use general input cfg for layout
layout = ft_prepare_layout(cfg);

%try replacing basic circle outline with extended outline when many
%channels
if numChan>128
    try
        cfg=[];
        cfg.layout='GSN-HydroCel-257.mat';

        layout_for_outline=ft_prepare_layout(cfg);

        layout.outline=layout_for_outline.outline;
    end
end

cfg=[];
cfg.xlim = [1 1];
cfg.zlim = [0 100]; %always scale between 0-100 for consistency
cfg.layout=layout;
cfg.style = 'straight'; %no contours
cfg.interplimits='electrode';
cfg.shading='flat'; %note: interp gives (visual) artifacts

%markers
cfg.marker='on';
cfg.markersymbol='.';

if numChan>=64
    cfg.markersize=4;
else
    cfg.markersize=8;
end

%bad channel markers
if ~isempty(bad_chan_inds) && strcmp(plot_style,'full')
    cfg.highlight ='on';
    cfg.highlightchannel=chanLabel(bad_chan_inds);
    cfg.highlightsymbol='.';
    cfg.highlightsize=8;
    cfg.highlightcolor='r';%artifact_grid_plot_cols(4,:);

end

%comment
cfg.comment='no';

%colorbar
cfg.colorbar = 'east';
cfg.colorbartext=[ grid_legend_labels{1} ' (%)']; % (dependent on plot_style)


%data
%trick fieldtrip by pretending we have chan x freq data
topo_data.powspctrm = 100*mean(artifact_data_topo,2); % (dependent on plot_style)
topo_data.freq = 1;
topo_data.label = chanLabel; % {1 x N}
topo_data.dimord = 'chan_freq';


ft_topoplotER(cfg,topo_data);

colormap(ax_topo,colormap_topo);

%adjust subplot
ax_topo.Position(2)=0.07;

%adjust colorbar
c=ax_topo.Colorbar;
c.Position(1)=0.9;
c.Position(4)=0.1;
set(c,'FontSize',fsize_legend,'Ticks',0:50:100)

%adjust colorbar label
l=c.Label;
set(l,'Position',[1.5 50 0],'Rotation',270)

%%
%titles
sgtitle(fh,strrep(fig_title,'_','\_'),'Fontsize',12,'FontWeight','bold')
fh.Name=fig_title;

% do the general cleanup and bookkeeping at the end of the function
% ft_postamble debug
% ft_postamble trackconfig
% ft_postamble previous data
% data = datanew
% ft_postamble provenance data
% ft_postamble history    data
% ft_postamble savevar    data

fprintf([functionname ' function finished\n']);
toc(ttic)
memtoc(mtic)
