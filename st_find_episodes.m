function episode_table=st_find_episodes(cfg,scoring)

epochs=scoring.epochs;%column orientation
nEpochs=length(epochs);


%set defaults


cfg.episodelabel=ft_getopt(cfg, 'episodelabel', 'none');
cfg.mincontinuouslength=ft_getopt(cfg,'mincontinuouslength',0.5);
cfg.interruptdetails              = ft_getopt(cfg, 'interruptdetails', {});
cfg.minepisodelength  = ft_getopt(cfg, 'minepisodelength', 5);

min_episode_length = floor(cfg.minepisodelength*60/scoring.epochlength);

%determine sleep onset/offset to restrict episodes
[onsetCandidateIndex, preOffsetCandidate, onsetepoch] = st_sleeponset(cfg,scoring);
if isempty(onsetepoch) % no sleep onset found
    onsetCandidateIndex = nEpochs;
end

if isempty(preOffsetCandidate)
    preOffsetCandidate = nEpochs;
end

%
consider_epochs=false(size(epochs));
consider_epochs(onsetCandidateIndex:preOffsetCandidate)=true;

%locate episode cores (continuous stretch meeting mincontinuouslength)
has_episode_label= ismember(epochs,cfg.episodestages) & consider_epochs;
has_episode_label_reps=getReps(has_episode_label);
episode_core=has_episode_label & has_episode_label_reps>=(cfg.mincontinuouslength*60/scoring.epochlength);
first_episode_core_start=find(episode_core==true,1,'first');
last_episode_core_end=find(episode_core==true,1,'last');

%set up episode gaps
episode_gaps=~episode_core;
if first_episode_core_start>1
    episode_gaps(1:first_episode_core_start-1)=false;
end
if last_episode_core_end<length(epochs)
    episode_gaps(last_episode_core_end+1:end)=false;
end
gap_change=diff([false episode_gaps false]);
episode_gap_onsets=find(gap_change==1);
episode_gap_offsets=find(gap_change==-1)-1;
episode_gap_inds=[episode_gap_onsets;episode_gap_offsets];
gap_allowed=true(1,size(episode_gap_inds,2));

%by default, set episode-terminating interruption as any mix of non-episode labels >0
maxinterruptdetails=cfg.interruptdetails;
if isempty(maxinterruptdetails)
    stages_present=unique(epochs);
    nonepisodestages=setdiff(stages_present,cfg.episodestages);

    maxinterruptdetails=[{nonepisodestages} num2cell(0)];
end

%loop across interruption labels (individual stages or sets of stages)

for interrupt_type_i=1:size(maxinterruptdetails,1)
    interruptlabel=maxinterruptdetails{interrupt_type_i,1};
    maxinterruptlength=maxinterruptdetails{interrupt_type_i,2};

    has_interrupt_label=ismember(epochs,interruptlabel) & episode_gaps;
    has_episode_interrupt_reps=getReps(has_interrupt_label);

    %locate sequences of interruption labels exceeding maximum length
    is_episode_ending_interrupt=has_interrupt_label & has_episode_interrupt_reps>(maxinterruptlength*60/scoring.epochlength);



    interrupt_starts=find(diff([false is_episode_ending_interrupt])==1);
    for interrupt_i=1:length(interrupt_starts)
        gap_allowed((episode_gap_inds(1,:)<=interrupt_starts(interrupt_i)) & (episode_gap_inds(2,:)>=interrupt_starts(interrupt_i)))=false;
    end

end

allowed_episode_gap_inds=episode_gap_inds(:,gap_allowed);
episode_filled=episode_core;

for fill_i=1:size(allowed_episode_gap_inds,2)
    episode_filled(allowed_episode_gap_inds(1,fill_i): allowed_episode_gap_inds(2,fill_i))=true;
end


fprintf('smoothing... %i epoch(s) altered\n',length(find(has_episode_label~=episode_filled)))

%start/end
% 
% %if first epoch is part of episode, already assign to starts here
% starts=[];
% if episode_filled(1) == 1
%     starts=1;
% end
% 
% %if last epoch is part of episode, add 0 for diff vector
% if episode_filled(end) == 1
%     episode_filled = [episode_filled 0];
% end

%determine episode starts/ends
onoff = diff([false episode_filled false]);
starts = find(onoff == 1);
ends   = find(onoff == -1)-1;

% if cfg.considersleeponset==true
%     %determine sleep onset/offset to restrict episodes
%     hasSleep = true;
%     [onsetCandidateIndex, preOffsetCandidate, onsetepoch] = st_sleeponset(cfg,scoring);
%     if isempty(onsetepoch) % no sleep onset found
%         hasSleep = false;
%         onsetCandidateIndex = nEpochs;
%     end
%
%     if isempty(preOffsetCandidate)
%         preOffsetCandidate = nEpochs;
%     end
%
%     %1) exclude episodes ending before sleep onset
%     idx = ends < onsetCandidateIndex;
%     starts(idx)=[];
%     ends(idx)=[];
%
%     %2) adjust episode starts occurring before sleep onset to sleep onset
%     starts(starts < onsetCandidateIndex) = onsetCandidateIndex;
%
% end

episode_lengths= ends-starts+1;

%3) exclude brief periods
idx=episode_lengths < min_episode_length;
starts(idx)=[];
ends(idx)=[];
episode_lengths(idx)=[];

%create output table
episode_table= table(...
    repmat({cfg.episodelabel},[numel(starts) 1]), [1:numel(starts)]',...
    starts(:), ends(:),episode_lengths(:),...
    'VariableNames',{...
    'episode_label','episode','startepoch', 'endepoch', 'durationepochs'}...
    );

    function Y=getReps(X)

        d = [true, diff(X) ~= 0, true];  % TRUE if values change
        n = diff(find(d));               % Number of repetitions
        Y = repelem(n, n);

    end

end