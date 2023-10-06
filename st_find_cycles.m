function [cycle_table,episode_table_all]=st_find_cycles(cfg,scoring)

cfg.allowsleeponsetbeforesleepopon= ft_getopt(cfg,'allowsleeponsetbeforesleepopon','no');

nEpochs=length(scoring.epochs);
%determine sleep onset/offset to restrict episodes
[onsetCandidateIndex, preOffsetCandidate, onsetepoch] = st_sleeponset(cfg,scoring);
if isempty(onsetepoch) % no sleep onset found
    onsetCandidateIndex = nEpochs;
end

if isempty(preOffsetCandidate)
    preOffsetCandidate = nEpochs;
end


%-------detect episodes-----
episode_table_all=table;




%W episodes
% cfg_episode=cfg;
% cfg_episode.episodelabel='epW';
% cfg_episode.episodestages={'W'};
% cfg_episode.mincontinuouslength=0.5;
% cfg_episode.interruptdetails={'R' 0;{'N1' 'N2' 'N3'} 2};
% cfg_episode.minepisodelength=5;
%
% episode_table=st_find_episodes(cfg_episode,scoring);
%episode_table_all=[episode_table_all;episode_table];

%W stable
cfg_episode=cfg;
cfg_episode.episodelabel='epW';
cfg_episode.episodestages={'W'};
cfg_episode.mincontinuouslength=2;
cfg_episode.interruptdetails={'R' 0;{'N1' 'N2' 'N3'} 1};
cfg_episode.minepisodelength=5;

episode_table=st_find_episodes(cfg_episode,scoring);
episode_table_all=[episode_table_all;episode_table];

%REM episodes
cfg_episode=cfg;
cfg_episode.episodelabel='epR';
cfg_episode.episodestages={'R'};
cfg_episode.mincontinuouslength=0.5;
cfg_episode.interruptdetails={'W' 2;{'W' 'N1' 'N2' 'N3'} 10};
cfg_episode.minepisodelength=0.5;

episode_table=st_find_episodes(cfg_episode,scoring);
episode_table_all=[episode_table_all;episode_table];

%NREM stab
cfg_episode=cfg;
cfg_episode.episodelabel='epNRst';
cfg_episode.episodestages={'N2','N3'};
cfg_episode.mincontinuouslength=2;
cfg_episode.interruptdetails={'W' 2;'N1' 2;'R' 0;{'W','N1'} 2};
cfg_episode.minepisodelength=20;

episode_table=st_find_episodes(cfg_episode,scoring);
episode_table_all=[episode_table_all;episode_table];

%NREM episodes
% cfg_episode=cfg;
% cfg_episode.episodelabel='epNR';
% cfg_episode.episodestages={'N1','N2','N3'};
% cfg_episode.mincontinuouslength=2;
% cfg_episode.interruptdetails={'W' 5;'R' 0};
% cfg_episode.minepisodelength=10;
%
% episode_table=st_find_episodes(cfg_episode,scoring);
%episode_table_all=[episode_table_all;episode_table];



%NREM frag
% cfg_episode=cfg;
% cfg_episode.episodelabel='epNRfr';
% % cfg_episode.episodestages={'N1','N2'};
% % cfg_episode.mincontinuouslength=0.5;
% %          cfg_episode.interruptdetails={'W' 5;'N3' 1;'R' 0;{'W','N3','R'} 5};
%
% cfg_episode.episodestages={'N1'};
% cfg_episode.mincontinuouslength=2;
% cfg_episode.interruptdetails={'W' 2;'N2' 2;'N3' 2;'R' 0;{'W','N2','N3'} 2};
% cfg_episode.minepisodelength=5;
%
% episode_table=st_find_episodes(cfg_episode,scoring);
%episode_table_all=[episode_table_all;episode_table];

%organize
episode_table_all=sortrows(episode_table_all,{'startepoch'});


% consider_epochs=false(size(scoring.epochs));
% episode_epochs=consider_epochs;
%
% consider_epochs(onsetCandidateIndex:preOffsetCandidate)=true;

%NREM frag
if ~isempty(episode_table_all)
    newEpLabel={'epNRtr'};
    epCount=1;
    newEps=table;
    %for each major episode, designate preceding gaps as epNRtr
    for ep_i=1:size(episode_table_all,1) +1

        epStart=[];
        if ep_i==1
            if onsetCandidateIndex<episode_table_all{ep_i,'startepoch'}
                epStart=onsetCandidateIndex;
                epEnd=episode_table_all{ep_i,'startepoch'}-1;

            end
        elseif ep_i==size(episode_table_all,1) +1
            if preOffsetCandidate>episode_table_all{ep_i-1,'endepoch'}
                epStart=episode_table_all{ep_i-1,'endepoch'}+1;
                epEnd=preOffsetCandidate;

            end

        else
            if episode_table_all{ep_i,'startepoch'}>episode_table_all{ep_i-1,'endepoch'}+1

                epStart=episode_table_all{ep_i-1,'endepoch'}+1;
                epEnd=episode_table_all{ep_i,'startepoch'}-1;

            end
        end


        if ~isempty(epStart)
            epDur=epEnd-epStart+1;
            newEp=table(newEpLabel,epCount,epStart,epEnd,epDur,...
                'VariableNames',{...
                'episode_label','episode','startepoch', 'endepoch', 'durationepochs'});

            newEps=[newEps;newEp];

            epCount=epCount+1;

        end




       

    end
     %add new episodes
        episode_table_all=[episode_table_all;newEps];
    %organize
    episode_table_all=sortrows(episode_table_all,{'startepoch'});


    %----from episodes to cycles---
    preWake_episodes=find(strcmp(episode_table_all.episode_label,'epW'))-1;
    R_episodes=find(strcmp(episode_table_all.episode_label,'epR'));
    cycle_end_WR_inds=union(preWake_episodes,R_episodes);

    lastEp=episode_table_all{end,'episode_label'};
    if ~ismember(lastEp,{'epW','epR'})
        cycle_end_WR_inds=union(cycle_end_WR_inds,size(episode_table_all,1)); %also include final episode
    end

    cycle_end_NR_inds=[];
    for end_WR_i=1:length(cycle_end_WR_inds)

        if end_WR_i==1
            startRowInd=1;
        else
            startRowInd=cycle_end_WR_inds(end_WR_i-1)+1;
        end
        endRowInd=cycle_end_WR_inds(end_WR_i);
        ind_local=find(strcmp(episode_table_all{startRowInd:endRowInd,'episode_label'},'epNRst'));
        if length(ind_local)>1
            cycle_end_NR_inds=[cycle_end_NR_inds; ind_local(1:end-1)+startRowInd-1];
        end

    end

    cycle_end_all_inds=union(cycle_end_WR_inds,cycle_end_NR_inds);

    %---now find coresponding cycle starts
    cycle_start_all_inds=[];
    for end_i=1:length(cycle_end_all_inds)

        if end_i==1
            cycle_start_ind=1;
        else
            cycle_start_ind=cycle_end_all_inds(end_i-1)+1;
            if strcmp(episode_table_all{cycle_start_ind,'episode_label'},'epW')
                cycle_start_ind=cycle_start_ind+1;
            end
        end
        cycle_start_all_inds=[cycle_start_all_inds;cycle_start_ind];
    end


    %build cycle table
    cycle_start_ends_inds=[cycle_start_all_inds(:) cycle_end_all_inds(:)];


    cycleNumber=[1:size(cycle_start_ends_inds,1)]';
    cycle_starts=episode_table_all{cycle_start_ends_inds(:,1),'startepoch'};
    cycle_ends=episode_table_all{cycle_start_ends_inds(:,2),'endepoch'};
    cycle_duration=cycle_ends-cycle_starts+1;

    cycle_composition=cellfun(@(X) strjoin(episode_table_all{X(1):X(2),'episode_label'},'-'),mat2cell(cycle_start_ends_inds,ones(size(cycle_start_ends_inds,1),1)),'UniformOutput',false);
    cycle_complete=contains(cycle_composition,'epNRst') & contains(cycle_composition,'epR');

    minCycleLength=30;
    cycle_complete=cycle_duration>=(minCycleLength*60/scoring.epochlength);

    cycle_table=table(cycleNumber,cycle_starts,cycle_ends,cycle_duration,cycle_composition,cycle_complete,...
        'VariableNames',{'cycle','startepoch','endepoch','durationepochs','episodes','complete'});
else
    cycle_table=table;

end

% switch episode_definition
%     case 'default' %15 min episode, 5 min W/R interruption allowed
%         cfg_episode.minepisodelength=15;
%         cfg_episode.maxinterruptlengthminutes=5;
%         cfg_episode.smoothexcludestages={};
%     case 'smart_interruptions' %15 min episode, 5 min W (not R) interruption allowed
%         cfg_episode.minepisodelength=15;
%         cfg_episode.maxinterruptlengthminutes=5;
%         cfg_episode.smoothexcludestages={'R'};
%     case 'sensitive' %15 min episode, 15 min W (not R) interruption allowed
%         cfg_episode.mincontinuouslength=10;
%          cfg_episode.interruptdetails={'W' 5;'R' 0};
%         cfg_episode.minepisodelength=15;
%
% end
%
% %find raw episodes
% episode_table=st_find_episodes(cfg_episode,scoring);
%
% %assign NREM labels to vector
% for ep_i=1:size(episode_table,1)
%     episode_labels(episode_table{ep_i,'startepoch'}:episode_table{ep_i,'endepoch'})=episode_table{1,'episode_label'};
% end
%
%
% episode_table_all=[episode_table_all;episode_table];
%
% %REM episodes
% cfg_episode=cfg;
% cfg_episode.episodelabel='R_ep';
% cfg_episode.episodestages={'R'};
% switch episode_definition
%     case 'default' %5 min episode, 0 min interruption allowed
%         cfg_episode.minepisodelength=5;
%         cfg_episode.smoothepisodeminutes=0;
%         cfg_episode.smoothexcludestages={};
%     case 'smart_interruptions' %5 min episode, 5 min W/NR interruption allowed
%         cfg_episode.minepisodelength=5;
%         cfg_episode.smoothepisodeminutes=5;
%         cfg_episode.smoothexcludestages={};
%     case 'sensitive' % 0.5 min episode, 15 min W/NR interruption allowed
%         cfg_episode.mincontinuouslength=0.5;
%         cfg_episode.interruptdetails={'W' 5 ;{'N1' 'N2' 'N3'} 5};
%         cfg_episode.minepisodelength=0.5;
% end
%
%
% %find raw episodes
% episode_table=st_find_episodes(cfg_episode,scoring);

% %assign REM labels to vector (MAY OVERWRITE NREM LABELS)
% for ep_i=1:size(episode_table,1)
%     episode_labels(episode_table{ep_i,'startepoch'}:episode_table{ep_i,'endepoch'})=episode_table{1,'episode_label'};
% end

% episode_table_all=[episode_table_all;episode_table];
% episode_table_all=sortrows(episode_table_all,{'startepoch'});

%%
% tmp_scoring=scoring;
% tmp_scoring.epochs=episode_labels';
%
% episode_table_all_2=table;
%
% %NREM episodes
% cfg_episode=cfg;
% cfg_episode.episodelabel='NR_ep';
% cfg_episode.episodestages={'NR_ep'};
% cfg_episode.minepisodelengthminutes=0.5;
% cfg_episode.smoothepisodeminutes=0;
% cfg_episode.smoothexcludestages={};
% cfg_episode.considersleeponset=false;
%
% episode_table=st_find_episodes(cfg_episode,tmp_scoring);
% episode_table_all_2=[episode_table_all_2;episode_table];
%
% %REM episodes
% cfg_episode=cfg;
% cfg_episode.episodelabel='R_ep';
% cfg_episode.episodestages={'R_ep'};
% cfg_episode.minepisodelengthminutes=0.5;
% cfg_episode.smoothepisodeminutes=0;
% cfg_episode.smoothexcludestages={};
% cfg_episode.considersleeponset=false;
%
% episode_table=st_find_episodes(cfg_episode,tmp_scoring);
% episode_table_all_2=[episode_table_all_2;episode_table];
%
% episode_table_all=sortrows(episode_table_all_2,{'startepoch'});


%create updated table
%

%[first_cycle_start,first_cycle_start_ind]=min(episode_table_all.startepoch)

% cycleStarts=[];
% cycleEnds=[];
%
% if ~isempty(episode_table_R) %end of R ep
%
%     %     first_cycle_end=episode_table_R{1,'endepoch'};
%     %
%     %     cycle_table=table(first_cycle_start,first_cycle_end,'VariableNames',{'startepoch','endepoch'});
%
%     for R_ep=1:size(episode_table_R)
%         if R_ep==1
%             cycleStarts = [cycleStarts;onsetCandidateIndex];
%         else
%             cycleStarts = [cycleStarts;  episode_table_R{R_ep-1,'endepoch'}+1]; %end of previous REM period (+1): could be W
%         end
%         cycleEnds   = [cycleEnds;  episode_table_R{R_ep,'endepoch'}];
%
%     end
%
% else %if no R ep, end of NR ep
%     first_cycle_end=episode_table{1,'endepoch'};
%
% end

% cycleLengths = (cycleEnds - cycleStarts + 1);
% startepochs = cycleStarts;
% endepochs = cycleEnds;
%
% maxCycle = max([numel(startepochs),numel(endepochs),numel(Rstarts),numel(Rends),numel(NRstarts),numel(NRends)]);
%
% startepochs((end+1):maxCycle) = NaN;
% endepochs((end+1):maxCycle) = NaN;
% Rstarts((end+1):maxCycle) = NaN;
% Rends((end+1):maxCycle) = NaN;
% NRstarts((end+1):maxCycle) = NaN;
% NRends((end+1):maxCycle) = NaN;
%
% for iCycle = 1:numel(startepochs)
%     if isnan(startepochs(iCycle))
%         startepochs(iCycle) = min([NRstarts(iCycle) Rstarts(iCycle)]);
%         endepochs(iCycle) = min([NRends(iCycle) Rends(iCycle)]);
%     end
% end
%
% startepochs = startepochs(:).'';
% endepochs = endepochs(:).'';
% Rstarts = Rstarts(:).'';
% Rends = Rends(:).'';
% NRstarts = NRstarts(:).'';
% NRends = NRends(:).'';



end

