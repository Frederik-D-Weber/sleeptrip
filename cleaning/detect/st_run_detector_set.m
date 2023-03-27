function cfg_artifacts=st_run_detector_set(cfg_detector_set,data)

% ST_RUN_DETECTOR_SET applies all detectors in cfg_detector_set (a detector set) to data,
% and returns details of individual artifacts in cfg_artifacts
%
% Use as:
%     cfg_artifacts=st_run_detector_set(cfg_detector_set,data)
%
% Required configuration parameters:
%     cfg.number      = numeric, number of detectors
%     cfg.label    = cell array of detector names
%     cfg.detectors = cell array of detector cfg structs
%
% Output:
%     cfg_artifacts = configuration containing channelwise artifact information
%
% See also ST_GET_DEFAULT_DETECTOR_SET

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

%TO DO:
%  consider adding mergeCloseIntervals and mergeOverlappingIntervals
%as subfunctions

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
%ft_preamble debug
%ft_preamble loadvar data
%ft_preamble provenance data
%ft_preamble trackconfig

% the ft_abort variable is set to true or false in ft_preamble_init
% if ft_abort
%   return
% end

%-----core function start---
%---input checks and defaults----
ft_checkconfig(cfg_detector_set,'required',{'number','label','detectors','elec'});

fprintf([functionname ' function initialized\n'])


%initialize cell array of eventTables
eventTables={};

[numChan, numSample]=size(data.trial{1});
srate=data.fsample;

%for every detector
detector_runtime=[];
for detector_i=1:cfg_detector_set.number

    tic_detector=tic;

    %get current detector
    cfg_detector=cfg_detector_set.detectors{detector_i};

    %verify existence of fields
    ft_checkconfig(cfg_detector,'required',{'label','st'});

    fprintf('detection artifacts of type: %s\n',cfg_detector.label);

    %----extract and check st cfg----
    cfg_st=cfg_detector.st;

    %verify existence of fields
    ft_checkconfig(cfg_st,'required',{'method'});

    %set required field defaults if not provided...
    cfg_st.minduration=ft_getopt(cfg_st,'minduration',0);
    cfg_st.maxduration=ft_getopt(cfg_st,'maxduration',inf);
    cfg_st.paddingduration=ft_getopt(cfg_st,'paddingduration',0);
    cfg_st.mergeduration=ft_getopt(cfg_st,'mergeduration',0);

    %...because we always need to determine these:
    minDurationSample=round(cfg_st.minduration*srate);
    maxDurationSample=round(cfg_st.maxduration*srate);
    padSample=round(cfg_st.paddingduration*srate);
    mergeSample=round(cfg_st.mergeduration*srate);



    % -------extract FieldTrip cfg and perform processing if requested
    if isfield(cfg_detector,'ft')
        cfg_ft=cfg_detector.ft;
        temp_data=ft_preprocessing(cfg_ft,data);
    else
        temp_data=data;
    end

    %get as simple matrix
    temp_data=temp_data.trial{1};

    %------perform SleepTrip processing
    if strcmp(cfg_st.method,'threshold')

        %verify existence of required threshold fields
        ft_checkconfig(cfg_st,'required',{'thresholddirection','thresholdvalue'});

        %Note the order in which signal processing is performed:
        %signal gradient if requested
        cfg_st.diff=ft_getopt(cfg_st,'diff','no');
        if strcmp(cfg_st.diff,'yes')
            temp_data=diff(temp_data, 1, 2);
        end

        %rectify if requested
        cfg_st.abs=ft_getopt(cfg_st,'abs','no');
        if strcmp(cfg_st.abs,'yes')
            temp_data=abs(temp_data);
        end

        %zscore if requested
        cfg_st.zscore=ft_getopt(cfg_st,'zscore','no');
        if strcmp(cfg_st.zscore,'yes')
            temp_data=zscore(temp_data,0,2); %across columns
        end

        %find samples meeting threshold
        if strcmp(cfg_st.thresholddirection,'above')
            dat_bad=temp_data>cfg_st.thresholdvalue;
        elseif strcmp(cfg_st.thresholddirection,'below')
            dat_bad=temp_data<cfg_st.thresholdvalue;
        elseif strcmp(cfg_st.thresholddirection,'between')
            dat_bad=temp_data>cfg_st.thresholdvalue(1) & temp_data<cfg_st.thresholdvalue(2);
        else
            ft_error('unknown option for thresholddirection')
        end

        %because of diff operation above, add first sample back as false (=
        %not artifact)
        if strcmp(cfg_st.diff,'yes')
            dat_bad=[false(numChan,1) dat_bad];
        end

        %get all start, end, and sample-duration values for each channel (using subfunction)
        start_end_sample=cellfun(@startEndSample, mat2cell(dat_bad,ones(numChan,1),numSample),'uni',false);


    elseif strcmp(cfg_st.method,'neighbours') %correlation with neighbors, below/above R and above proportion of channels

        %verify existence of required fields
        ft_checkconfig(cfg_st,'required',{'elec','neighbours','metric','thresholddirection','thresholdvalue'});

        %segment length, get default  and process
        cfg_st.segmentlength=ft_getopt(cfg_st,'segmentlength',30); %default 30 s segmentlength
        segmentlength_sample=round(cfg_st.segmentlength*srate);
        numSegment  = floor(numSample/segmentlength_sample); %excludes final, incomplete segment (handled later)

        %get connectivity matrix an number of connections per channel (using neighbours and channel labels)
        cfg_st.channel=cfg_st.elec.label;
        connectivity=channelconnectivity(cfg_st);
        numConnected=sum(connectivity,2);

        %initialize matrices
        badchannel = false(numChan,numSegment);
        chanChanSeg=zeros(numChan,numChan,numSegment); %for plotting/debugging only



        for segment_i=1:numSegment

            %get segment start/ends
            segment_start_sample=(segment_i-1)*segmentlength_sample+1;
            segment_end_sample=segment_i*segmentlength_sample;

            %if final segment, set end sample to end of data (=include incomplete final segment)
            if segment_i==numSegment
                segment_end_sample=numSample;
            end

            %segment data
            dat_seg=temp_data(:,segment_start_sample:segment_end_sample);


            % functions for calculating all-to-all metrics
            if strcmp(cfg_st.metric,'correlation')
                %channel correlation
                chanChanMetric=corr(dat_seg');
                chanChanMetric(~connectivity)=nan;

            elseif strcmp(cfg_st.metric,'maxabsdiff')
                %channel max magnitude of difference wave
                chanChanMetric=maxAbsDiff(dat_seg');
                chanChanMetric(~connectivity)=nan;

            else
                ft_error('unknown option for metric')

            end


            %store for plotting/debugging
            chanChanSeg(:,:,segment_i)=chanChanMetric;

            %find channels meeting threshold
            if strcmp(cfg_st.thresholddirection,'above')
                %dat_bad=temp_data>cfg_st.thresholdvalue;
                chanchan_bad=chanChanMetric>cfg_st.thresholdvalue;
            elseif strcmp(cfg_st.thresholddirection,'below')
                chanchan_bad=chanChanMetric<cfg_st.thresholdvalue;
            else
                ft_error('unknown option for thresholddirection')
            end

            %label segment bad if proportion of channels exceeds channelthreshold
            badchannel(:,segment_i) = (sum(chanchan_bad,2)./numConnected)>cfg_st.channelthreshold;


        end

        %find idices of bad segments
        badTrial_inds=cellfun(@(X) find(X), mat2cell(badchannel,ones(numChan,1),numSegment),'uni',false);

        %get all start, end, and sample-duration values
        start_end_sample=cellfun(@(X) [(X'-1)*segmentlength_sample+1 X'*segmentlength_sample repmat(segmentlength_sample,[size(X,2) 1])],...
            badTrial_inds,'uni',false);
    else
        ft_error('unknown option for method')
    end


    %---process start_end_sample times---

    %select the events exceeding minDurationSample
    start_end_sample=cellfun(@(X) X(X(:,3)>minDurationSample & X(:,3)<maxDurationSample,:),...
        start_end_sample,'uni',false);

    %pad with requested length padSample (not retaining original duration)
    start_end_sample=cellfun(@(X) [X(:,1)-padSample X(:,2)+padSample],...
        start_end_sample,'uni',false);

    %merge overlapping intervals
    start_end_sample=cellfun(@mergeOverlappingIntervals,start_end_sample,'uni',false);

    %merge nearby events based on requested definition
    %Note: mergeSample+1 so intervals immediately adjacent also get merged
    %Note: adding "intervals" of begin and end sample to allow extension to
    %data start/end if it falls within merge interval
    start_end_sample=cellfun(@(X) mergeCloseIntervals([1 1; X; numSample numSample],mergeSample+1),...
        start_end_sample,'uni',false);

    %remove intervals with identical start/end
    start_end_sample=cellfun(@(X) X(X(:,1)~=X(:,2),:),start_end_sample,'uni',false);


    %-------build table----
    %add channel label
    start_end_sample=cellfun(@(X,Y) [num2cell(X) repmat({Y},[size(X,1) 1])],...
        start_end_sample,data.label,'uni',false);

    %convert to table, merge channels
    eventTable=cell2table(vertcat(start_end_sample{:}),'VariableNames',{'start' 'stop' 'channel'});

    %adjust events outside data range
    eventTable{eventTable{:,'start'}<1,'start'}=1;
    eventTable{eventTable{:,'stop'}>numSample,'stop'}=numSample;

    %convert to time in s
    %     eventTable.start=(eventTable.start-1)./srate;
    %     eventTable.stop=(eventTable.stop-1)./srate;

    %convert to time in seconds (respecting time vector in data)
    eventTable.start=(eventTable.start-1)./srate+data.time{1}(1);
    eventTable.stop=(eventTable.stop-1)./srate+data.time{1}(1);

    %add duration
    eventTable.duration=eventTable.stop-eventTable.start+1/srate;

    eventTable(:,'event')={cfg_detector.label};

    %reorder
    [~, varOrder] = ismember(eventTable.Properties.VariableNames, {'event','channel','start','stop','duration'});
    [~, resortOrder] = sort(varOrder);
    eventTable = eventTable(:,resortOrder);

    %add to cell array of eventTables
    eventTables{detector_i}=eventTable;

    detector_runtime(detector_i)=toc(tic_detector);
end

cfg_artifacts=[];
cfg_artifacts.continuous.label=cfg_detector_set.label;
cfg_artifacts.continuous.artifacts=eventTables;
cfg_artifacts.detector_set=cfg_detector_set;
cfg_artifacts.data=rmfield(data,{'trial'}); %add data struct (without actual data)
cfg_artifacts.elec=cfg_detector_set.elec;
cfg_artifacts.detector_runtime=detector_runtime;

% do the general cleanup and bookkeeping at the end of the function
%ft_postamble debug
%ft_postamble trackconfig
% ft_postamble previous data
% data = datanew
% ft_postamble provenance data
% ft_postamble history    data
% ft_postamble savevar    data

fprintf([functionname ' function finished\n']);
toc(ttic)
memtoc(mtic)

end

%--------------SUBFUNCTIONS---------

function chanChanMetric=maxAbsDiff(data)

%the actual function computing difference waves (expand into 3rd dim to prevent element-wise subtraction):
metricFunc=@(X) squeeze(max(abs(X-reshape(X,size(X,1),1,size(X,2)))));

chanChanMetric=metricFunc(data);

end

%for a boolean vector, find start/end indices and duration of TRUE stretches
function start_end_sample=startEndSample(X)
startSamples=find(diff([0 X])==1); %precede with 0 so data starting with artifact is caught
endSamples= find(diff([X 0])==-1); %follow with 0 so data ending with artifact is caught

start_end_sample=[startSamples' endSamples' (endSamples-startSamples)'+1]; %+1, so single sample has duration of 1
end


