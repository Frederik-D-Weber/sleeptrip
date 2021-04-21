function [f_k,f_sek,f_p,f_z,f_ci,f_kj,f_sekj,f_zkj,f_pkj,...
          consensusModusVoteMatrix,consensusCountMatrix,cross_comparisons,chi2_cross,p_cross, ...
          ctbls,scoring_num_ref_labels,scoring_num_to_ref_labels] = mulitple_kappa(alignmentMatrix,alphalevel,referenceOption)
% Copyright (C) 2019-, Frederik D. Weber 
%               with code from Giuseppe Cardillo
%               with code from  Cardillo G. (2007) Fleiss'es kappa: compute the Fleiss'es kappa for multiple raters.   
%               adapted from Giuseppe Cardillo giuseppe.cardillo-edta@poste.it
%               https://www.mathworks.com/matlabcentral/fileexchange/15426-fleiss
%               https://github.com/dnafinder/Fleiss
%               Cardillo G. (2007) Cohen's kappa: compute the Cohen's kappa ratio on a square matrix.   
%               http://www.mathworks.com/matlabcentral/fileexchange/15365
%               https://github.com/dnafinder/Cohen
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

% alignmentMatrix = [0 1 1 1 2   2 3 3 6 1 NaN;...
%                    0 1 1 1 1   2 2 3 6 1 2  ;...
%                    0 1 1 1 NaN 2 2 3 5 1 NaN]
categories = unique(alignmentMatrix);
categories = categories(~isnan(categories));
nCategories = numel(categories);
%referenceOption = 'first'; %either 'first' or 'consenus'

nObservations = size(alignmentMatrix,2);
nSamples = size(alignmentMatrix,1);


%%% Fleiss'es kappa

consensusCountMatrix = zeros(nCategories,nObservations);
consensusModusVoteMatrix = zeros(1,nObservations);

for iColumn = 1:nObservations
    for iCat = 1:nCategories
        consensusCountMatrix(iCat,iColumn) = sum(alignmentMatrix(:,iColumn) == categories(iCat));
    end
    consensusModusVoteMatrix(iColumn) = categories(find(consensusCountMatrix(:,iColumn) == max(consensusCountMatrix(:,iColumn)),1,'first'));
end

if nCategories > 1
[f_k,f_sek,f_p,f_z,f_ci,f_kj,f_sekj,f_zkj,f_pkj] = fleiss(consensusCountMatrix',alphalevel);
else
  f_k = 1;
  f_sek = 0;
  f_p = 0;
  f_z = NaN;
  f_ci = repmat(nan,1,nSamples);
  f_kj =  repmat(1,1,nSamples);
  f_sekj =  repmat(1,1,nSamples);
  f_zkj =  repmat(nan,1,nSamples);
  f_pkj =  repmat(0,1,nSamples);
end

switch referenceOption
    case 'first'
        reference_index = 1;
        reference_label = '1';
        reference = alignmentMatrix(reference_index,:);
        startCompIndex = 2;
    case 'consensus'
        reference_index = -1;
        reference_label = 'consensus';
        reference = consensusModusVoteMatrix;
        startCompIndex = 1;
    case 'all'
        reference_index = -1;
    otherwise
        ft_error('referenceOption called %s unknown.',referenceOption)
end




cross_comparisons = {};
chi2_cross = [];
p_cross = [];
nComparison = 0;
scoring_num_ref_labels = {};
scoring_num_to_ref_labels = {};
switch referenceOption
    case {'first' 'consensus'}
        for iComp = startCompIndex:nSamples
            nComparison = nComparison+1;
            
            test = alignmentMatrix(iComp,:);
            [tbl,chi2_cross_temp,p_cross_temp,labels] = crosstab(test,reference);
            row_lab = labels(:,1);
            row_lab_index = cellfun(@(s) ~isempty(s),row_lab,'UniformOutput',false);
            row_lab_index = [row_lab_index{:}];
            row_lab = row_lab(row_lab_index);
            
            col_lab = labels(:,2);
            col_lab_index = cellfun(@(s) ~isempty(s),col_lab,'UniformOutput',false);
            col_lab_index = [col_lab_index{:}];
            col_lab = col_lab(col_lab_index);
            
            row_lab = cellfun(@(s) strrep(s,'-1','unscored'), row_lab,'UniformOutput',false);
            col_lab = cellfun(@(s) strrep(s,'-1','unscored'), col_lab,'UniformOutput',false);
            
%             row_lab = cellfun(@(s) strrep(s,'-2','unassigned'), row_lab,'UniformOutput',false);
%             col_lab = cellfun(@(s) strrep(s,'-2','unassigned'), col_lab,'UniformOutput',false);
            
            cross_comparisons{nComparison} = array2table(tbl,'RowNames',cellfun(@(s) ['comp_' s], row_lab,'UniformOutput',false),'VariableNames',cellfun(@(s) ['ref_' s], col_lab,'UniformOutput',false));
            chi2_cross(nComparison) = chi2_cross_temp;
            p_cross(nComparison) = p_cross_temp;
            scoring_num_ref_labels = cat(1,scoring_num_ref_labels,{reference_label});
            scoring_num_to_ref_labels = cat(1,scoring_num_to_ref_labels,{num2str(iComp)});
        end
    case 'all'
        for iComp = 1:(nSamples-1)
            reference_label = num2str(iComp);
            for iComp2 = (iComp+1):nSamples
                nComparison = nComparison+1;
                
                reference = alignmentMatrix(iComp,:);
                test = alignmentMatrix(iComp2,:);
                [tbl,chi2_cross_temp,p_cross_temp,labels] = crosstab(test,reference);
                row_lab = labels(:,1);
                row_lab_index = cellfun(@(s) ~isempty(s),row_lab,'UniformOutput',false);
                row_lab_index = [row_lab_index{:}];
                row_lab = row_lab(row_lab_index);
                
                col_lab = labels(:,2);
                col_lab_index = cellfun(@(s) ~isempty(s),col_lab,'UniformOutput',false);
                col_lab_index = [col_lab_index{:}];
                col_lab = col_lab(col_lab_index);
                
                row_lab = cellfun(@(s) strrep(s,'-1','unscored'), row_lab,'UniformOutput',false);
                col_lab = cellfun(@(s) strrep(s,'-1','unscored'), col_lab,'UniformOutput',false);
                
%                 row_lab = cellfun(@(s) strrep(s,'-2','unassigned'), row_lab,'UniformOutput',false);
%                 col_lab = cellfun(@(s) strrep(s,'-2','unassigned'), col_lab,'UniformOutput',false);
                
                cross_comparisons{nComparison} = array2table(tbl,'RowNames',cellfun(@(s) ['comp_' s], row_lab,'UniformOutput',false),'VariableNames',cellfun(@(s) ['ref_' s], col_lab,'UniformOutput',false));
                chi2_cross(nComparison) = chi2_cross_temp;
                p_cross(nComparison) = p_cross_temp;
                scoring_num_ref_labels = cat(1,scoring_num_ref_labels,{reference_label});
                scoring_num_to_ref_labels = cat(1,scoring_num_to_ref_labels,{num2str(iComp2)});
            end
        end
end



%%% Cohen's kappa

ctbls = [];
nComparison = 0;

%for iComp = 1:nComparison
%     alignmentMatrix_pair = cat(1,alignmentMatrix(iComp,:),reference);
%
%     categories_pair = unique(alignmentMatrix_pair);
%     categories_pair = categories_pair(~isnan(categories_pair));
%     nCategories_pair = numel(categories_pair);
%
%     nObservations_pair = size(alignmentMatrix_pair,2);
%     nSamples_pair = size(alignmentMatrix_pair,1);
%
%
%     consensusCountMatrix_pair = zeros(nCategories_pair,nObservations_pair);
%
%     for iColumn = 1:nObservations_pair
%         for iCat = 1:nCategories_pair
%             consensusCountMatrix_pair(iCat,iColumn) = sum(alignmentMatrix_pair(:,iColumn) == categories_pair(iCat));
%         end
%     end

switch referenceOption
    case {'first' 'consensus'}
        for iComp = startCompIndex:nSamples
            nComparison = nComparison+1;
            
            test = alignmentMatrix(iComp,:);
            confusionmat_pair = confusionmat(reference,test);
            
            [c_k,c_sek,c_ci,c_kmax,c_kratiomax,c_po,c_pe,c_po_minus_pe,c_one_minus_pe,c_alpha,c_vari,c_z,c_p] = cohenskappa(confusionmat_pair, alphalevel);
            
            pair = cat(1,reference,test);
            agreement = sum(all(bsxfun(@eq,pair,pair(1,:))))/size(pair,2);
            
            tbl_first = table(scoring_num_ref_labels(nComparison),scoring_num_to_ref_labels(nComparison),...
                'VariableNames',{'scoring_num_ref','scoring_num_comp'});
            
            
            tbl_second = array2table([agreement,c_k,c_sek,c_p,c_z,c_ci(1),c_ci(2),c_kmax,c_kratiomax,c_po,c_pe,c_po_minus_pe,c_one_minus_pe,c_alpha,c_vari],...
                'VariableNames',{'agreement','Cohens_kappa','Cohens_kappa_SEM','Cohens_kappa_pvalue','Cohens_kappa_zvalue','Cohens_kappa_CI_lower','Cohens_kappa_CI_higher','Cohens_kappa_max_possible','Cohens_kappa_ratio_of_max_possible','agreement_observed','agreement_expected','agreement_percentage_due_to_true_concordance','agreement_percentage_residual_not_random','alpha_level','Cohens_kappa_variance'});
            
            tbl = cat(2,tbl_first,tbl_second);
            
            if isempty(ctbls)
                ctbls = tbl;
            else
                ctbls = cat(1,ctbls,tbl);
            end
        end
    case 'all'
        for iComp = 1:(nSamples-1)
            reference_label = num2str(iComp);
            for iComp2 = (iComp+1):nSamples
                nComparison = nComparison+1;
                
                reference = alignmentMatrix(iComp,:);
                test = alignmentMatrix(iComp2,:);
                
                confusionmat_pair = confusionmat(reference,test);
                
                [c_k,c_sek,c_ci,c_kmax,c_kratiomax,c_po,c_pe,c_po_minus_pe,c_one_minus_pe,c_alpha,c_vari,c_z,c_p] = cohenskappa(confusionmat_pair, alphalevel);
                
                pair = cat(1,reference,test);
                agreement = sum(all(bsxfun(@eq,pair,pair(1,:))))/size(pair,2);
                
                
                tbl_first = table(scoring_num_ref_labels(nComparison),scoring_num_to_ref_labels(nComparison),...
                    'VariableNames',{'scoring_num_ref','scoring_num_comp'});
                
                tbl_second = array2table([agreement,c_k,c_sek,c_p,c_z,c_ci(1),c_ci(2),c_kmax,c_kratiomax,c_po,c_pe,c_po_minus_pe,c_one_minus_pe,c_alpha,c_vari],...
                    'VariableNames',{'agreement','Cohens_kappa','Cohens_kappa_SEM','Cohens_kappa_pvalue','Cohens_kappa_zvalue','Cohens_kappa_CI_lower','Cohens_kappa_CI_higher','Cohens_kappa_max_possible','Cohens_kappa_ratio_of_max_possible','agreement_observed','agreement_expected','agreement_percentage_due_to_true_concordance','agreement_percentage_residual_not_random','alpha_level','Cohens_kappa_variance'});
                tbl = cat(2,tbl_first,tbl_second);
                
                if isempty(ctbls)
                    ctbls = tbl;
                else
                    ctbls = cat(1,ctbls,tbl);
                end
            end
        end
end



%end

end



function [k,sek,p,z,ci,kj,sekj,zkj,pkj] = fleiss(varargin)
% cite: Cardillo G. (2007) Fleiss'es kappa: compute the Fleiss'es kappa for multiple raters.   
% adapted from Giuseppe Cardillo giuseppe.cardillo-edta@poste.it
% https://www.mathworks.com/matlabcentral/fileexchange/15426-fleiss
% https://github.com/dnafinder/Fleiss
% FLEISS: compute the Fleiss'es kappa
% Fleiss'es kappa is a generalisation of Scott's pi statistic, a
% statistical measure of inter-rater reliability. It is also related to
% Cohen's kappa statistic. Whereas Scott's pi and Cohen's kappa work for
% only two raters, Fleiss'es kappa works for any number of raters giving
% categorical ratings (see nominal data), to a fixed number of items. It
% can be interpreted as expressing the extent to which the observed amount
% of agreement among raters exceeds what would be expected if all raters
% made their ratings completely randomly. Agreement can be thought of as
% follows, if a fixed number of people assign numerical ratings to a number
% of items then the kappa will give a measure for how consistent the
% ratings are. The scoring range is between 0 and 1.
%
% Syntax: 	fleiss(X,alpha)
%
%     Inputs:
%           X - square data matrix
%           ALPHA - significance level (default = 0.05)
%     Outputs:
%           - kappa value for the j-th category (kj)
%           - kj standard error
%           - z of kj
%           - p-value
%           - Fleiss'es kappa
%           - kappa standard error
%           - kappa confidence interval
%           - k benchmarks by Landis and Koch
%           - z test
%
%      Example:
%An example of the use of Fleiss'es kappa may be the following: Consider
%fourteen psychiatrists are asked to look at ten patients. Each
%psychiatrist gives one of possibly five diagnoses to each patient. The
%Fleiss' kappa can be computed from this matrix to show
%the degree of agreement between the psychiatrists above the level of
%agreement expected by chance.
%x =
%     0     0     0     0    14
%     0     2     6     4     2
%     0     0     3     5     6
%     0     3     9     2     0
%     2     2     8     1     1
%     7     7     0     0     0
%     3     2     6     3     0
%     2     5     3     2     2
%     6     5     2     1     0
%     0     2     2     3     7
%
% So there are 10 rows (1 for each patient) and 5 columns (1 for each
% diagnosis). Each cell represents the number of raters who
% assigned the i-th subject to the j-th category
% x=[0 0 0 0 14; 0 2 6 4 2; 0 0 3 5 6; 0 3 9 2 0; 2 2 8 1 1 ; 7 7 0 0 0;...
% 3 2 6 3 0; 2 5 3 2 2; 6 5 2 1 0; 0 2 2 3 7];
%
%  Calling on Matlab the function: fleiss(x);
%
%           Answer is:
%
% kj:       0.2013    0.0797    0.1716    0.0304    0.5077
%
% s.e.:     0.0331
%
% z:        6.0719    2.4034    5.1764    0.9165   15.3141
%
% p:        0.0000    0.0162    0.0000    0.3594         0
%
% ------------------------------------------------------------
% Fleiss'es (overall) kappa = 0.2099
% kappa error = 0.0170
% kappa C.I. (95%) = 0.1767 	 0.2432
% Fair agreement
% z = 12.3743 	 p = 0.0000
% Reject null hypotesis: observed agreement is not accidental
%
%           Created by Giuseppe Cardillo
%           giuseppe.cardillo-edta@poste.it
%
% To cite this file, this would be an appropriate format:
% Cardillo G. (2007) Fleisses kappa: compute the Fleiss'es kappa for multiple raters.
% http://www.mathworks.com/matlabcentral/fileexchange/15426

args=cell(varargin);
nu=numel(args);
if isempty(nu)
    error('Warning: Matrix of data is missed...')
elseif nu>2
    error('Warning: Max two input data are required')
end
default.values = {[],0.05};
default.values(1:nu) = args;
[x alpha] = deal(default.values{:});

if isvector(x)
    error('Warning: x must be a matrix, not a vector.');
end
if ~all(isfinite(x(:))) || ~all(isnumeric(x(:))) || isempty(x)
    error('Warning: X data matrix values must be numeric and finite')
end
if ~isequal(x(:),round(x(:)))
    error('Warning: X data matrix values must be whole numbers')
end
n=size(x,1); %subjects
%chech if the raters are the same for each rows
r=sum(x,2);
if any(r-max(r))
    error('The raters are not the same for each rows')
end
if nu==2 %if necessary check alpha
    if ~isscalar(alpha) || ~isnumeric(alpha) || ~isfinite(alpha) || isempty(alpha)
        error('Warning: it is required a numeric, finite and scalar ALPHA value.');
    end
    if alpha <= 0 || alpha >= 1 %check if alpha is between 0 and 1
        error('Warning: ALPHA must be comprised between 0 and 1.')
    end
end
clear args default nu

m=sum(x(1,:)); %raters
a=n*m;
pj=(sum(x)./(a)); %overall proportion of ratings in category j
b=pj.*(1-pj);
c=a*(m-1);
d=sum(b);
kj=1-(sum((x.*(m-x)))./(c.*b)); %the value of kappa for the j-th category
sekj=realsqrt(2/c); %kj standar error
zkj=kj./sekj;
pkj=(1-0.5*erfc(-abs(zkj)/realsqrt(2)))*2;
k=sum(b.*kj)/d; %Fleiss'es (overall) kappa
sek=realsqrt(2*(d^2-sum(b.*(1-2.*pj))))/sum(b.*realsqrt(c)); %kappa standard error
ci=k+([-1 1].*(abs(0.5*erfc(-alpha/2/realsqrt(2)))*sek)); %k confidence interval
z=k/sek; %normalized kappa
p=(1-0.5*erfc(-abs(z)/realsqrt(2)))*2;

% %display results
% fprintf('kj:   '); disp(kj)
% fprintf('s.e.: '); disp(sekj)
% fprintf('z:    '); disp(zkj)
% fprintf('p:    '); disp(pkj)
% disp(repmat('-',1,60))
% fprintf('Fleiss''es (overall) kappa = %0.4f\n',k)
% fprintf('kappa error = %0.4f\n',sek)
% fprintf('kappa C.I. (%d%%) = %0.4f \t %0.4f\n',(1-alpha)*100,ci)
% if k<0
%     disp('Poor agreement')
% elseif k>=0 && k<=0.2
%     disp('Slight agreement')
% elseif k>0.2 && k<=0.4
%     disp('Fair agreement')
% elseif k>0.4 && k<=0.6
%     disp('Moderate agreement')
% elseif k>0.6 && k<=0.8
%     disp('Substantial agreement')
% elseif k>0.8 && k<=1
%     disp('Perfect agreement')
% end
% fprintf('z = %0.4f \t p = %0.4f\n',z,p)
% if p<0.05
%     disp('Reject null hypotesis: observed agreement is not accidental')
% else
%     disp('Accept null hypotesis: observed agreement is accidental')
% end
end


function [k,sek,ci,kmax,kratiomax,po,pe,po_minus_pe,one_minus_pe,alpha,vari,z,p] = cohenskappa(x,alpha)
% KAPPA: This function computes the Cohen's kappa coefficient.
% Cohen's kappa coefficient is a statistical measure of inter-rater
% reliability. It is generally thought to be a more robust measure than
% simple percent agreement calculation since k takes into account the
% agreement occurring by chance.
% Kappa provides a measure of the degree to which two judges, A and B,
% concur in their respective sortings of N items into k mutually exclusive
% categories. A 'judge' in this context can be an individual human being, a
% set of individuals who sort the N items collectively, or some non-human
% agency, such as a computer program or diagnostic test, that performs a
% sorting on the basis of specified criteria.
% The original and simplest version of kappa is the unweighted kappa
% coefficient introduced by J. Cohen in 1960. When the categories are
% merely nominal, Cohen's simple unweighted coefficient is the only form of
% kappa that can meaningfully be used. 
%
% Syntax: 	kappa(X,ALPHA)
%      
%     Inputs:
%           X - square data matrix
%           ALPHA - default=0.05.
%
%     Outputs:
%           - Observed agreement percentage
%           - Random agreement percentage
%           - Agreement percentage due to true concordance
%           - Residual not random agreement percentage
%           - Cohen's kappa 
%           - kappa error
%           - kappa confidence interval
%           - Maximum possible kappa
%           - k observed as proportion of maximum possible
%           - k benchmarks by Landis and Koch 
%           - z test results
%
%      Example: 
%
%           x=[88 14 18; 10 40 10; 2 6 12];
%
%           Calling on Matlab the function: kappa(x)
%
%           Answer is:
%
% UNWEIGHTED COHEN'S KAPPA
% --------------------------------------------------------------------------------
% Observed agreement (po) = 0.7000
% Random agreement (pe) = 0.4100
% Agreement due to true concordance (po-pe) = 0.2900
% Residual not random agreement (1-pe) = 0.5900
% Cohen's kappa = 0.4915
% kappa error = 0.0549
% kappa C.I. (alpha = 0.0500) = 0.3839     0.5992
% Maximum possible kappa, given the observed marginal frequencies = 0.8305
% k observed as proportion of maximum possible = 0.5918
% Moderate agreement
% Variance = 0.0031     z (k/sqrt(var)) = 8.8347    p = 0.0000
% Reject null hypotesis: observed agreement is not accidental
%
%           Created by Giuseppe Cardillo
%           giuseppe.cardillo-edta@poste.it
%
% To cite this file, this would be an appropriate format:
% Cardillo G. (2007) Cohen's kappa: compute the Cohen's kappa ratio on a square matrix.   
% http://www.mathworks.com/matlabcentral/fileexchange/15365
m=size(x,1);

f=diag(ones(1,m)); %unweighted
n=sum(x(:)); %Sum of Matrix elements
x=x./n; %proportion
r=sum(x,2); %rows sum
s=sum(x); %columns sum
Ex=r*s; %expected proportion for random agree
pom=sum(min([r';s])); %maximum proportion observable
po=sum(sum(x.*f)); %proportion observed
pe=sum(sum(Ex.*f)); %proportion expected

k=(po-pe)/(1-pe); %Cohen's kappa
kmax=(pom-pe)/(1-pe); %maximum possible kappa, given the observed marginal frequencies
kratiomax=k/kmax; %observed as proportion of maximum possible
sek=sqrt((po*(1-po))/(n*(1-pe)^2)); %kappa standard error for confidence interval
ci=k+([-1 1].*(abs(-realsqrt(2)*erfcinv(alpha))*sek)); %k confidence interval
wbari=r'*f;
wbarj=s*f;
wbar=repmat(wbari',1,m)+repmat(wbarj,m,1);
a=Ex.*((f-wbar).^2);
vari=(sum(a(:))-pe^2)/(n*((1-pe)^2)); %variance
z=k/sqrt(vari); %normalized kappa
p=(1-0.5*erfc(-abs(z)/realsqrt(2)))*2;

po_minus_pe = (po-pe);
one_minus_pe = (1-pe);

%         %display results
%         fprintf('Observed agreement (po) = %0.4f\n',po)
%         fprintf('Random agreement (pe) = %0.4f\n',pe)
%         fprintf('Agreement due to true concordance (po-pe) = %0.4f\n',po-pe)
%         fprintf('Residual not random agreement (1-pe) = %0.4f\n',1-pe)
%         fprintf('Cohen''s kappa = %0.4f\n',k)
%         fprintf('kappa error = %0.4f\n',sek)
%         fprintf('kappa C.I. (alpha = %0.4f) = %0.4f     %0.4f\n',alpha,ci)
%         fprintf('Maximum possible kappa, given the observed marginal frequencies = %0.4f\n',kmax)
%         fprintf('k observed as proportion of maximum possible = %0.4f\n',kratiomax)
%         if k<0
%             disp('Poor agreement')
%         elseif k>=0 && k<=0.2
%             disp('Slight agreement')
%         elseif k>=0.21 && k<=0.4
%             disp('Fair agreement')
%         elseif k>=0.41 && k<=0.6
%             disp('Moderate agreement')
%         elseif k>=0.61 && k<=0.8
%             disp('Substantial agreement')
%         elseif k>=0.81 && k<=1
%             disp('Perfect agreement')
%         end
%         fprintf('Variance = %0.4f     z (k/sqrt(var)) = %0.4f    p = %0.4f\n',vari,z,p)
%         if p<alpha
%             disp('Reject null hypotesis: observed agreement is not accidental')
%         else
%             disp('Accept null hypotesis: observed agreement is accidental')
%         end
end