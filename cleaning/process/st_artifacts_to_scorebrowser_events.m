function eventTable=st_artifacts_to_scorebrowser_events(cfg_artifacts)

ft_checkconfig(cfg_artifacts,'required',{'artifacts','scorebrowser_events','data'});



eventTable={};

%locate requested event types in raw_events
if isfield(cfg_artifacts.artifacts,'raw_events')
    eventTypesAvailable=fieldnames(cfg_artifacts.artifacts.raw_events);
    includeEventTypes=intersect(cfg_artifacts.scorebrowser_events,eventTypesAvailable);

    for evType_i=1:length(includeEventTypes)
        events=cfg_artifacts.artifacts.raw_events.(includeEventTypes{evType_i}).events;


        if ~isempty(events)
            events=eventTable_expand_all(events,cfg_artifacts.data.label);
        end

        eventTable{end+1}=events;
    end
end

%locate requested event types in grid_events
if isfield(cfg_artifacts.artifacts,'grid_events')
    eventTypesAvailable=fieldnames(cfg_artifacts.artifacts.grid_events);
    includeEventTypes=intersect(cfg_artifacts.scorebrowser_events,eventTypesAvailable);

    for evType_i=1:length(includeEventTypes)
        events=cfg_artifacts.artifacts.grid_events.(includeEventTypes{evType_i}).events;

        if ~isempty(events)
            events=eventTable_expand_all(events,cfg_artifacts.data.label);
        end

        eventTable{end+1}=events;
    end
end

%concatenate
eventTable=vertcat(eventTable{:});