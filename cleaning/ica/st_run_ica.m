function [comp,compLabelProbsTable]=st_run_ica(cfg,data)


cfg.preclean=ft_getopt(cfg, 'preclean', 'yes'); %clean by default
cfg.stages_ica_train=ft_getopt(cfg, 'stages_ica_train', {'W','N1','R'});
cfg.ica_opts=ft_getopt(cfg, 'ica_opts', {});


if istrue(cfg.preclean)

    %get standard detector set if nothing supplied
    if ~isfield(cfg,'detector_set')
        cfg_tmp=[];
        cfg_tmp.elec=cfg.elec;
        cfg_tmp.fsample=data.fsample;
        cfg_tmp.include={'highamp'};
        cfg_detector_set=st_get_default_detector_set(cfg_tmp);
    end

    %use user-supplied detector set, otherwise default for ICA
    cfg.detector_set  = ft_getopt(cfg, 'detector_set', cfg_detector_set);

    %detect artifacts
    cfg_artifacts=st_run_detector_set(cfg.detector_set,data);

    %process the artifacts (including by stage)
    cfg_artifacts.scoring=cfg.scoring;
    cfg_artifacts.segmentrejectthresh=0.5; %exclude segment if half of channels are artifactual
    cfg_artifacts=st_process_detector_results(cfg_artifacts);

    %
    scoring_for_selection=cfg_artifacts.scoring_artifact_level;
else
    scoring_for_selection=cfg.scoring;
    %chan_prop_bad=mean(cfg_artifacts.artifacts.grid.artifact_grid_segment_expanded,2);
end

%select desired stages
cfg_select=[];
cfg_select.scoring=scoring_for_selection;
cfg_select.stages=cfg.stages_ica_train;
cfg_select.usescoringexclusion='yes';
data_ica=st_select_data(cfg_select, data);

%make continuous
data_ica=st_trial_to_continuous([],data_ica);


%extract data
dat=data_ica.trial{1};
data_minutes=round(size(dat,2)/data.fsample)/60;

fprintf('%.1f min of data available for ICA...\n', data_minutes)
if data_minutes<60
    ft_warning('less than an hour of data available for ICA...\n')
end



%run ICA (with options from cfg.ica_opts)
[weights, sphere] = runica(dat, cfg.ica_opts{:});

% calculate mixing/unmixing matrices
unmixing = weights * sphere; %comp x chan
mixing = pinv(unmixing); %chan x comp


%----initialize EEGLAB struct---
EEG=eeg_emptyset;
EEG.srate=data_ica.fsample;
EEG.data=dat; %concatenated data
EEG.chanlocs=elec2chanlocs(cfg.elec); %custom func
EEG.nbchan= length(data.label);
EEG.trials=1;
EEG.pnts=size(dat,2);

%ICA (icaact not needed: eeg_checkset will remove)
EEG.icachansind=1:length(data.label);
EEG.icaweights=weights;
EEG.icasphere=sphere;
EEG.icawinv=mixing;

EEG=eeg_checkset(EEG);

%call IClabel
EEG=iclabel(EEG); %will rereference to average, recalculate component time courses


%get classification tables
ic_classification=EEG.etc.ic_classification.ICLabel;
[m,maxInd]=max(ic_classification.classifications,[],2);
compLabelTable=table(ic_classification.classes(maxInd)',m,'VariableNames',{'classLabel','prob'});
compLabelProbsTable=array2table(ic_classification.classifications,'VariableNames',ic_classification.classes);



%---Fieldtrip---

cfg_comp = [];
cfg_comp.demean    = 'no';           % This has to be explicitly stated, as the default is to demean.
cfg_comp.unmixing  = unmixing;  % Supply the matrix necessary to 'unmix' the channel-series data into components
cfg_comp.topolabel = data.label(:); % Supply the original channel label information

%calculate component time courses for entire data
comp = ft_componentanalysis(cfg_comp, data);

%adjust component names
for k = 1:size(comp.topo,2)
    comp.label{k,1} = sprintf('%03d_%s_%.2f', k,compLabelTable{k,'classLabel'}{1}(1:3),compLabelTable{k,'prob'});
end


%--remove  components


%     cfg=[];
%     cfg.component=badCompInds;
%     data=ft_rejectcomponent(cfg, comp);
%
%     %simplify output comp structure (can be restored with ft_componentanalysis)
%     comp=rmfield(comp,{'time','trial'});
%
%     badcomp=[];
%     badcomp.inds=badCompInds;
%     badcomp.probs=compLabelProbsTable;

