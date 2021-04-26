function [times, nEvents, property] = eventSmoother(eventsSeconds,eventsProperties,windowInSeconds,timeStepInSeconds,startTime,endTime)

times = startTime:timeStepInSeconds:endTime;
nEvents = zeros(1,numel(times));
property = zeros(1,numel(times));

for iTimeStep = 1:numel(times);
    timeStep = times(iTimeStep);
    leftBoundary  = timeStep-windowInSeconds/2;
    rightBoundary = timeStep+windowInSeconds/2;
    eventIndex = find((eventsSeconds >= leftBoundary) & (eventsSeconds < rightBoundary));
    nEvents(iTimeStep)  = numel(eventIndex);
    property(iTimeStep) = mean(eventsProperties(eventIndex));
end
end