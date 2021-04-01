function [tempfinaldestinationpaths tempfinaldestinationfolder] = st_irodsget(cfg)

% ST_IRODSGET gets files (via iget command) from an irods server evironment into a local
% temporary folder for local useage
%
% Use as
%   [tempfinaldestinationpaths tempfinaldestinationfolder] = st_irodsget(cfg)
%
% Required configuration parameters are:
%   cfg.irodspaths = a string (single path) or a cellstr (multiple paths) with the paths as they apear on
%                    the irods file system 
%
% Optional configuration parameters are:
%   cfg.tempfolderbase  = string with a folder base, only change if you want 
%                         to manually clean up later without defaults of ST_IRODSCLEANUP 
%                         (default = '~/st_irodsget')
%
%
% See also ST_IRODSCLEANUP, ST_WAIT_EXISTS

% Copyright (C) 2021-, Frederik D. Weber
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
ft_preamble init
ft_preamble debug
%ft_preamble loadvar data
%ft_preamble provenance data
ft_preamble trackconfig

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  return
end

% set defaults
cfg.tempfolderbase = ft_getopt(cfg, 'tempfolderbase', ['~/' functionname ]); 

fprintf([functionname ' function initialized\n']);

if ~isfield(cfg, 'irodspaths')
    ft_error('cfg.irodspath needs to be defined as a full path in the irods environment.')
end

if isempty(cfg.irodspaths)
    ft_error('cfg.irodspath needs to be not empty.')
end

if ~iscellstr(cfg.irodspaths)
    cfg.irodspaths = {cfg.irodspaths};
end


%create uuid
tempuuid =  java.util.UUID.randomUUID;
tempuuid = tempuuid.toString;
uuid = tempuuid.toCharArray';


%create datetime string
dt = now;
timestampstring = [datestr(dt,'yyyy-mm-dd-HH-MM-SS-FFF')];

%here the files will be parked 
tempuniquesubfolder = [timestampstring '-' uuid];

tempfinaldestinationfolder = [cfg.tempfolderbase '/' tempuniquesubfolder];
tempfinaldestinationfolder_restartfile = [cfg.tempfolderbase '/' tempuniquesubfolder 'iget_restartfile.txt'];

[status_mkdir, message_mkdir, messageid_mkdir] = mkdir(tempfinaldestinationfolder);

if status_mkdir ~= 1
    ft_error('failed to make temporary folder to get files: %s',message_mkdir)
end

tempfinaldestinationpaths = {};
for iirodspaths = 1:numel(cfg.irodspaths)
irods_pathrequested = icfg.irodspaths{iirodspaths};

% remove a trailing path separator from irods as iget anyway ignores this!
if strcmp(irods_pathrequested(end),'/')
    irods_pathrequested = irods_pathrequested(1:(end-1));
end

command = ['iget -K -r -P -N 10 -X ' tempfinaldestinationfolder_restartfile ' --retries 10 ' irods_pathrequested ' ' tempfinaldestinationfolder];

[status_system, cmdout_system] = system(command);

if status_system ~= 0
    ft_error('failed to execute %d of %d iget commands: %s',iirodspaths, numel(cfg.irodspaths), command)
end

[pathstr, name, ext] = fileparts(irods_pathrequested);
finaldestinationfile = [tempfinaldestinationfolder '/' name ext];
tempfinaldestinationpaths = cat(1, tempfinaldestinationpaths, finaldestinationfile);

end

if numel(tempfinaldestinationpaths) == 1
    tempfinaldestinationpaths = tempfinaldestinationpaths{1};
end
%res:
%tempfinaldestinationpaths
%tempfinaldestinationfolder

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble trackconfig
% ft_postamble previous data
% data = datanew
% ft_postamble provenance data
% ft_postamble history    data
% ft_postamble savevar    data

fprintf([functionname ' function finished\n']);
toc(ttic)
memtoc(mtic)
end
