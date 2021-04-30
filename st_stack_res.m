function [table] = st_stack_res(varargin)

% ST_STACK_RES stack result structures by converting the results into a long format table and
% then appending them (stack-append)
%
% Use as
%   [table] = st_stack_res(res,...)
%
%
% See also ST_APPEND_RES

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
memtic
functionname = getfunctionname();
fprintf([functionname ' function started\n']);

% check if some are empty


empty_res = [];
iE = 0;
for iArg = 1:nargin
    if isempty(varargin{iArg})
        iE = iE + 1;
        empty_res(iE) = iArg;
        ft_warning('Result number %d is empty and is thus ignored. Continous numbering is also ignoring this result',iArg)
    end
end

if ~isempty(empty_res)
    varargin(empty_res) = [];
end




nRes = numel(varargin);

resIDs = cell(nargin,1);


r = varargin{1};
o = r.ori;
t = r.type;

allAppended = true;
anyAppended = false;

for iRes = 1:nRes
    r = varargin{iRes};
    

    if ~isempty(r.table)
    numericVars = varfun(@isnumeric,r.table,'output','uniform');
    nonNumericVarNames = r.table.Properties.VariableNames(find(~numericVars));
    if isempty(nonNumericVarNames)
        descriptor = repmat({'?'},size(r.table,1),1);
    else
        descriptor = strcat(nonNumericVarNames{1},'_',r.table.(nonNumericVarNames{1}));
        if numel(nonNumericVarNames) > 1
            for iNNVN = 2:(numel(nonNumericVarNames))
                descriptor = strcat(descriptor,'_',nonNumericVarNames{iNNVN},'_',r.table.(nonNumericVarNames{iNNVN}));
            end
        end
    end
    wasAppended = false;
    if isfield(r,'appended')
        wasAppended = r.appended;
    end
    
    if wasAppended
        idx_resnum = find(ismember(r.table.Properties.VariableNames, {'resnum'}),'first',1);
        if ~isempty(idx_resnum)
            rescol = r.table(:,idx_resnum);
            if isnumeric(rescol)
                rescol = num2str(table2array(rescol));
            end
            	descriptor = strcat('resnum','_',rescol,'_',descriptor);
        end
    end
        
    rt = cat(2,table(repmat(iRes,numel(descriptor),1),repmat({r.ori},numel(descriptor),1),repmat({r.type},numel(descriptor),1),descriptor,'VariableNames',{'resstacknum','res_ori','res_type','descriptor'}),r.table(:,find(numericVars)));
    rt = stack(rt,5:size(rt,2),...
                  'IndexVariableName','property',...
                  'NewDataVariableName','value');
    else
    
    tempIDnames = cat(2,{'resstacknum'},{'res_ori'},{'res_type'},{'descriptor'},{'property'},{'value'});
    tempvarnames = [tempIDnames];
    rt = cell2table(cell(0,numel(tempvarnames)), 'VariableNames', tempvarnames);
    end
    
    anyAppended = anyAppended || wasAppended;
    allAppended = allAppended && wasAppended;
end
    

if anyAppended
   ft_warning('some results have been appended before, will put the resnum column in the descriptor as well for the results affected')
end

fprintf([functionname ' function finished\n']);
toc(ttic)
memtoc
end