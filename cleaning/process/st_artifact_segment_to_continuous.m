function cfg_artifacts=st_artifact_segment_to_continuous(cfg_artifacts)

%---input checks and defaults----
ft_checkconfig(cfg_artifacts,'required',{'data'});
ft_checkconfig(cfg_artifacts.artifacts,'required',{'grid'});

cfg_grid=cfg_artifacts.artifacts.grid;


segmentLength=cfg_grid.segment_length;
chanLabels=cfg_artifacts.channel;


sampleinfo=cfg_artifacts.data.sampleinfo;
srate=cfg_artifacts.data.fsample;
segmentLengthSample=round(segmentLength*srate);

grid_names=fieldnames(cfg_grid)';
grid_names=grid_names(contains(grid_names,'grid'));

%initialize artifacts struct
artifacts=struct;
for grid_var_i=1:length(grid_names)

    eventLabel=grid_names{grid_var_i};

    eventLabelSimple=['grid_' erase(eventLabel,{'grid_','artifact_','_grid'})]; %for later use in event tables and artifacts.grid_events
    segments=cfg_grid.(eventLabel);

    %only consider 2D grids
    if ~isequal(size(segments),[cfg_grid.channel_number cfg_grid.segment_number])
        continue
    end

    %for grids below, only use first row of grids
    if ismember(eventLabel,{'reject_grid','accept_grid'})
        segments=segments(1,:);
        useChanLabels={'all'};
    else
        useChanLabels=chanLabels;
    end

    %
    dataDims=size(segments);


    %get each start/end sample
    start_end_sample=cellfun(@(X) [(find(X)-1)*segmentLengthSample+1; find(X)*segmentLengthSample]',...
        mat2cell(segments,ones(dataDims(1),1),dataDims(2)),'uni',false);


    %merge nearby intervals ()
    start_end_sample=cellfun(@(X) mergeCloseIntervals(X,1),start_end_sample,'uni',false);

    %add channel label
    start_end_sample=cellfun(@(X,Y) [num2cell(X) repmat({Y},[size(X,1) 1])],start_end_sample,useChanLabels,'uni',false);

    %convert to table, merge channels
    eventTable=cell2table(vertcat(start_end_sample{:}),'VariableNames',{'start' 'stop' 'channel'});


    %adjust events outside data range
    eventTable{eventTable{:,'start'}<1,'start'}=1;
    eventTable{eventTable{:,'stop'}>sampleinfo(2),'stop'}=sampleinfo(2);

    %convert to time in s
    eventTable.start=(eventTable.start-1)./srate;
    eventTable.stop=(eventTable.stop-1)./srate;


    %add duration
    eventTable.duration=eventTable.stop-eventTable.start;

    eventTable(:,'event')={eventLabelSimple};

    %reorder
    [~, varOrder] = ismember(eventTable.Properties.VariableNames, {'event','channel','start','stop','duration'});
    [~, resortOrder] = sort(varOrder);
    eventTable = eventTable(:,resortOrder);


    %collect grid event artifact info
    artifacts.(eventLabelSimple).label=eventLabelSimple;
    artifacts.(eventLabelSimple).events=eventTable;


end

cfg_artifacts.artifacts.grid_events=artifacts;
