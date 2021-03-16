function [res_comp_fleiss_stat res_comp_cohen_stat res_contingency res_contingency_excluded  scoring_consensus contingency_tables contingency_excluded_tables] = st_scoringcomp(cfg, varargin)
% ST_SCORINGCOMP compares several scorings against each other with measures of
% accuracy and Fleiss's as well as Cohen's  kappa as well as a consensus (in case of ties the
% scorings given first will prevail.) and the contingency tables with
% respect to the reference. If the reference is the first scoring, then
% the first item in the cel of contingency_tables output is empty.
% also see http://www.tqmp.org/RegularArticles/vol08-1/p023/ for a paper on
% the subject.
%
% Use as
%   [res_comp_fleiss_stat res_comp_cohen_stat res_contingency res_contingency_excluded  scoring_consensus contingency_tables contingency_excluded_tables] = st_scoringcomp(cfg, scoring1, scoring2, ...)
%
% Optional configuration parameters are:
%
%   cfg.sleeponsetdef   = string, sleep onset either 'N1' or 'N1_NR' or 'N1_XR' or
%                         'NR' or 'N2R' or 'XR' or 'AASM' or 'X2R' or
%                         'N2' or 'N3' or 'SWS' or 'S4' or 'R',
%                         see below for details (default = 'N1_XR') see
%                         ST_SLEEPONSET for details
%   cfg.align           = if scorings should be aligned to an epoch
%                         'sleeponset' or 'lighsoff' or 'no'
%                         or a single epoch number (default = 'no')
%                         see ST_CUTSCORING cfg.start parameter for details
%   cfg.stat_alpha      = value for statistical alpha level (default is 0.05)
%   cfg.reference       = which of the scorings should be chosen as reference 
%                         either 'first' or 'consensus' (default = 'first')
%   cfg.agreementthres  = number for the minimal agreement threshold in the range [0 1] 
%                         to decide if the consensus scoring in an epoch is reached. 
%                         0.5 means half of the raters needs to agree to
%                         reach consensus.
%                         (default = 0.5)
%   cfg.agreementthres_excl = number for the minimal agreement threshold in the range [0 1] 
%                         to decide if the consensus scoring if an epoch is excluded is reached. 
%                         0.5 means half of the raters needs to agree to
%                         reach consensus to exclude this epoch.
%                         (default = 0.5)
%
% See also ST_READ_SCORING, ST_SCORINGCONVERT, ST_SLEEPONSET, ST_CUTSCORING

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

% set defaults
%cfg.channel  = ft_getopt(cfg, 'channel', 'all', 1);
cfg.align          = ft_getopt(cfg, 'align', 'no');
cfg.stat_alpha     = ft_getopt(cfg, 'stat_alpha', 0.05);
cfg.reference      = ft_getopt(cfg, 'reference', 'first');
cfg.agreementthres = ft_getopt(cfg, 'agreementthres', 0.5);
cfg.agreementthres_excl = ft_getopt(cfg, 'agreementthres_excl', 0.5);


%TODO: consider scoring offsets that might differ,

scorings = varargin;

if ~iscell(scorings)
    scorings = {scorings};
end

nScorings = numel(scorings);
if nScorings < 2
    ft_error('At least two scorings expected but only given %d.',nScorings)
end


for iScoring = 2:nScorings
    scoring = scorings{iScoring};
    if ~(strcmp(scoring.standard,scorings{1}.standard) && ...
            (scoring.epochlength == scorings{1}.epochlength) && ...
            (scoring.dataoffset == scorings{1}.dataoffset))
        ft_error('Mixed scoring standards, epoch lengths or data offsets not supported yet!\n To unify scoring standards consider using ST_SCORINGCONVERT')
    end
end


numbersNorm = [];
excludedNorm = [];

for iScoring = 1:nScorings
    scoring = scorings{iScoring};
    
    
    if ~strcmp(cfg.align,'no')
        cfg.start = cfg.align;
        scorings_temp = st_cutscoring(cfg, scoring);
        scorings{iScoring} = scorings_temp{1};
    end
    
    cfg_scc = [];
    cfg_scc.to = 'number';
    scoring_numbers = st_scoringconvert(cfg_scc,scoring);
    scoring.numbers = cellfun(@str2num,scoring_numbers.epochs,'UniformOutput',1);
    
    if iScoring > 1
        nEpochs = size(numbersNorm,2);
        if numel(scoring.epochs) > nEpochs
            missingEpochs =  numel(scoring.epochs) - nEpochs;
            numbersNorm   =  [numbersNorm  NaN(size(numbersNorm,1),  missingEpochs)];
            excludedNorm  =  [excludedNorm NaN(size(excludedNorm,1), missingEpochs)];
        end
        
        if numel(scoring.epochs) < nEpochs
            missingEpochs = nEpochs - numel(scoring.epochs);
            scoring.numbers = [scoring.numbers, NaN(1,missingEpochs)];
            scoring.excluded = [scoring.excluded, NaN(1,missingEpochs)];
        end
    end
    numbersNorm(iScoring,:) = scoring.numbers;
    excludedNorm(iScoring,:) = scoring.excluded;
    
end


numbersNorm(isnan(numbersNorm)) = -1;

[k_norm,sek_norm,p_norm,z_norm,ci_norm,kj_norm,sekj_norm,zkj_norm,pkj_norm,...
    consensusModusVoteMatrix_norm,consensusModusCountMatrix_norm,cross_comparisons_norm,chi2_cross_norm,p_cross_norm,...
    ctbls] = ...
    multiple_kappa(numbersNorm,cfg.stat_alpha,cfg.reference);

% rename the labels to the ones that fit with the original scoring
for iCrossComp = 1:numel(cross_comparisons_norm)
    %cross_comparisons_norm{2}
    
    ccn = cross_comparisons_norm{iCrossComp};
    
    if ~isempty(ccn)
        
        scoremap = [];
        scoremap.unknown   = 'ref_unknown';
        scoremap.labelnew  = {'0', '1', '2', '3', '4', '5', '-1', '8'};
        scoring_dummy = scorings{1};
        scoring_dummy.excluded = [];
        scoring_dummy.standard = 'custom';
        cfg_scc = [];
        cfg_scc.to = 'custom';
        
        scoring_label_conti = {'ref_0', 'ref_1', 'ref_2', 'ref_3', 'ref_4', 'ref_5', 'ref_unscored', 'ref_8'};
        scoring_dummy.epochs = ccn.Properties.VariableNames;
        scoring_dummy.label = scoring_label_conti;
        scoremap.labelold  = scoring_label_conti;
        cfg_scc.scoremap = scoremap;
        scoring_dummy_ref = st_scoringconvert(cfg_scc,scoring_dummy);
        scoring_dummy_ref.standard = 'number';
        
        scoring_label_conti = {'comp_0', 'comp_1', 'comp_2', 'comp_3', 'comp_4', 'comp_5', 'comp_unscored', 'comp_8'};
        scoring_dummy.epochs = ccn.Properties.RowNames;
        scoring_dummy.label = scoring_label_conti;
        scoremap.labelold  = scoring_label_conti;
        cfg_scc.scoremap = scoremap;
        
        scoring_dummy_comp = st_scoringconvert(cfg_scc,scoring_dummy);
        scoring_dummy_comp.standard = 'number';
        %scoring_dummy.label
        
        cfg_scc = [];
        %cfg_scc.scoremap = scoremap;
        cfg_scc.to = scorings{1}.standard;
        scoring_dummy_ref = st_scoringconvert(cfg_scc,scoring_dummy_ref);
        scoring_dummy_comp = st_scoringconvert(cfg_scc,scoring_dummy_comp);
        
        ref_labels = scoring_dummy_ref.epochs;
        comp_labels = scoring_dummy_comp.epochs;
        
        ref_labels = cellfun(@(s) ['ref_' s], ref_labels,'UniformOutput',false);
        comp_labels = cellfun(@(s) ['comp_' s], comp_labels,'UniformOutput',false);
        
        ref_labels = cellfun(@(s) strrep(s,'?','unscored'), ref_labels,'UniformOutput',false);
        comp_labels = cellfun(@(s) strrep(s,'?','unscored'), comp_labels,'UniformOutput',false);
        
        ccn.Properties.VariableNames = ref_labels;
        ccn.Properties.RowNames = comp_labels;
        
        cross_comparisons_norm{iCrossComp} = ccn;
        
    end
end

excludedNorm(isnan(excludedNorm)) = -1;

[k_MA,sek_MA,p_MA,z_MA,ci_MA,kj_MA,sekj_MA,zkj_MA,pkj_MA,consensusModusVoteMatrix_MA,consensusModusCountMatrix_MA,cross_comparisons_MA,chi2_cross_MA,p_cross_MA...
          ctbls_MA] = ...
    multiple_kappa(excludedNorm,cfg.stat_alpha,cfg.reference);

agreement_epochs_single = all(bsxfun(@eq,numbersNorm,numbersNorm(1,:)));
agreement_excluded_single = all(bsxfun(@eq,excludedNorm,excludedNorm(1,:)));

agreement = sum(agreement_epochs_single)/size(numbersNorm,2);
agreement_excluded = sum(agreement_excluded_single)/size(excludedNorm,2);



temp_nRows = 1;%numel(kj_norm);
temp_stats = repmat([agreement,agreement_excluded,k_norm,k_MA,sek_norm,sek_MA,p_norm,p_MA,z_norm,z_MA,ci_norm(1),ci_norm(2),ci_MA(1),ci_MA(2)],temp_nRows,1);
temp_stats_table = array2table(temp_stats,'VariableNames',{'agreement','agreement_excluded','Fleiss_kappa','Fleiss_kappa_excluded','Fleiss_kappa_SEM','Fleiss_kappa_SEM_excluded','Fleiss_kappa_pvalue','Fleiss_kappa_excluded_pvalue','Fleiss_kappa_zvalue','Fleiss_kappa_excluded_zvalue','Fleiss_kappa_CI_lower','Fleiss_kappa_CI_higher','Fleiss_kappa_excluded_CI_lower','Fleiss_kappa_excluded_CI_higher'});

res_comp_fleiss_stat = [];
res_comp_fleiss_stat.ori = functionname;
res_comp_fleiss_stat.type = 'comp_fleiss_stats';
res_comp_fleiss_stat.cfg = cfg;
res_comp_fleiss_stat.table = temp_stats_table;


res_comp_cohen_stat = [];
res_comp_cohen_stat.ori = functionname;
res_comp_cohen_stat.type = 'comp_cohen_stats';
res_comp_cohen_stat.cfg = cfg;

if ~isempty(ctbls_MA)
    ctbls_MA.Properties.VariableNames = cellfun(@(s) [s '_excluded'], ctbls_MA.Properties.VariableNames,'UniformOutput',false);
end
res_comp_cohen_stat.table = cat(2,ctbls, ctbls_MA);

contingency_tables = cross_comparisons_norm;
contingency_excluded_tables = cross_comparisons_MA;



res_contingency = [];
res_contingency.ori = functionname;
res_contingency.type = 'contingency';
res_contingency.cfg = cfg;
res_contingency.table = concatCrossComp(cross_comparisons_norm,chi2_cross_norm,p_cross_norm);


res_contingency_excluded = [];
res_contingency_excluded.ori = functionname;
res_contingency_excluded.type = 'contingency_excluded';
res_contingency_excluded.cfg = cfg;
res_contingency_excluded.table = concatCrossComp(cross_comparisons_MA,chi2_cross_MA,p_cross_MA);

%make a new consensus scoring from ground up and exclude fields
%that are not crucial
scoringnew = [];
scoringnew.label = scorings{1}.label;
scoringnew.epochlength = scorings{1}.epochlength;
scoringnew.numbers = consensusModusVoteMatrix_norm;
scoringnew.epochs = cellfun(@num2str,num2cell(scoringnew.numbers),'UniformOutput',0);
idxNoConsensusReached = ~any((consensusModusCountMatrix_norm .* 1/nScorings) >= cfg.agreementthres);
scoringnew.numbers(idxNoConsensusReached) = -1;
scoringnew.epochs(idxNoConsensusReached) = {'-1'};
%scoring_consensus.prob = scorings_cycle_prop_adjusted;
scoringnew.excluded = logical(consensusModusVoteMatrix_MA);
idxNoConsensusReached = ~any((consensusModusCountMatrix_MA .* 1/nScorings) >= cfg.agreementthres_excl);
scoringnew.excluded(idxNoConsensusReached) = false;
scoringnew.standard = 'number';
scoringnew.agreement = max(consensusModusCountMatrix_norm .* 1/nScorings);
scoringnew.agreement_excluded = max(consensusModusCountMatrix_MA .* 1/nScorings);
scoringnew.full_agreement = logical(agreement_epochs_single);
scoringnew.full_agreement_excluded = logical(agreement_excluded_single);
scoringnew.dataoffset = scorings{1}.dataoffset;

cfg_scc = [];
cfg_scc.to = scoring.standard;
scoring_consensus = st_scoringconvert(cfg_scc, scoringnew);

fprintf([functionname ' function finished\n']);
toc(ttic)
memtoc(mtic)
end


function cc_all = concatCrossComp(cross_comparisons,chi2_cross,p_cross)
colnames = {};
for iCrossComp = 1:numel(cross_comparisons)
    colnames = cat(2,colnames,cross_comparisons{iCrossComp}.Properties.VariableNames);
end

colnames = unique(colnames);

for iCrossComp = 1:numel(cross_comparisons)
    tmptbl = cross_comparisons{iCrossComp};
    missing = find(~ismember(colnames,tmptbl.Properties.VariableNames));
    for iMissing = missing
        tmptbl.(colnames{iMissing}) = zeros(size(tmptbl,1),1);
    end
    cross_comparisons{iCrossComp} = tmptbl(:,colnames);
end

cc_all = [];
for iCrossComp = 1:numel(cross_comparisons)
    %iCrossComp = 2
    if isempty(cross_comparisons{iCrossComp})
        continue;
    end
    temp_cc_nrows = size(cross_comparisons{iCrossComp},1);
    
    temp_scoringnumber_table = array2table(repmat(iCrossComp,temp_cc_nrows,1));
    temp_stage_code_table = cell2table(cross_comparisons{iCrossComp}.Properties.RowNames);
    temp_chi2_table = array2table(repmat(chi2_cross(iCrossComp),temp_cc_nrows,1));
    temp_chi2_pvalue_table = array2table(repmat(p_cross(iCrossComp),temp_cc_nrows,1));
    
    
    temp_scoringnumber_table.Properties.VariableNames = {'scoring_num'};
    temp_stage_code_table.Properties.VariableNames = {'stage_code_comp_scoring'};
    temp_chi2_table.Properties.VariableNames = {'comp_chi2'};
    temp_chi2_pvalue_table.Properties.VariableNames = {'comp_chi2_pvalue'};
    
    
    temp_cc = [temp_scoringnumber_table temp_chi2_table temp_chi2_pvalue_table temp_stage_code_table cross_comparisons{iCrossComp}];
    
    temp_cc.Properties.RowNames =  cellstr(arrayfun(@(x) [num2str(iCrossComp) '_' x], num2str((1:size(temp_cc,1))'),'UniformOutput',false));
    
    if isempty(cc_all)
        cc_all = temp_cc;
    else
        cc_all = [cc_all ; temp_cc];
    end
end
end
