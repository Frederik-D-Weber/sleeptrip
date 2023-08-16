function eventTable=eventTable_expand_all(eventTable,labels)

events_all_inds=strcmp(eventTable.channel,'all');

eventTable_all=eventTable(events_all_inds,:);
eventTable(events_all_inds,:)=[];

numEventsAll=size(eventTable_all,1);
numLabels=length(labels);

eventTable_all_expanded=[];
for ev_i=1:numEventsAll
    temp=repmat(eventTable_all(ev_i,:),[numLabels 1]);
    temp.channel=labels;
    
    eventTable_all_expanded=[eventTable_all_expanded; temp];
end

eventTable=[eventTable; eventTable_all_expanded];