function undone = st_pipe_progress(id,process,done)

% ST_PIPE_PROGRESS keeps track of a the process (a string) with a mat file unique to the
% id (a string) and optionally sets if it is done (a logical). it returns
% if the process requested or set is undone;
%
% Use as
%   undone = st_pipe_progress(id,process) % to request/init the process status
%   undone = st_pipe_progress(id,process,done) % to set/overwrite the process status
%
% See also ST_WAIT_EXISTS

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

% these are used by the ft_preamble/ft_postamble function and scripts
ft_revision = '$Id$';
ft_nargin   = nargin;
ft_nargout  = nargout;

% do the general setup of the function
st_defaults
% ft_preamble init
% ft_preamble debug
% ft_preamble loadvar data
% ft_preamble provenance data
% ft_preamble trackconfig

% the ft_abort variable is set to true or false in ft_preamble_init
%if ft_abort
%  return
%end

% set defaults
%cfg.channel  = ft_getopt(cfg, 'channel', 'all', 1);
                                    
fprintf([functionname ' function initialized\n']);

undone = true;
doSave = false;
set_done = true;
if nargin < 3
    done = false;
    set_done = false;
end
filename = ['_' 'st_pipe' '-' id '.mat'];

if (exist(filename) == 2)
    load(st_wait_exists([],filename),'pipe')
else
    pipe = [];
    pipe.process = {};
    pipe.done = {};
end

idx_status = ismember(pipe.process,process);
if isempty(idx_status) || ~any(idx_status)
    pipe.process = cat(2,pipe.process,{process});
    pipe.done = cat(2,pipe.done,{done});
    undone = ~done;
    doSave = true;
else
    if (done ~= pipe.done{idx_status})
        if set_done
            pipe.done{idx_status} = done;
            doSave = true;
        end
    end
    undone = ~pipe.done{idx_status};
end

if doSave
    save(filename,'pipe');
end


% do the general cleanup and bookkeeping at the end of the function
% ft_postamble debug
% ft_postamble trackconfig
% ft_postamble previous data
% data = datanew
% ft_postamble provenance data
% ft_postamble history    data
% ft_postamble savevar    data

fprintf([functionname ' function finished\n']);
toc(ttic)
memtoc(mtic)
end