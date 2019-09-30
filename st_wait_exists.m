function [filename] = st_wait_exists(cfg, filename, varargin)

% ST_WAIT_EXISTS pauses program until filename exists
%
% Use as
%   [waited_seconds] = st_wait_load(cfg, filename, varargin)
%
%   cfg can be empty (cfg = [];)
%   varargin is passed to 'load' function
%
% Optional configuration parameters are:
%
%   cfg.timeout         = scalar, time in seconds to wait at most (default = Inf)
%   cfg.checkinterval   = scalar, interval in seconds to check (default = 1)
%   cfg.writerateapprox = scalar, assumed write speed in Bytes per second of a modified file (default = 1000000)
%   cfg.minmodwait      = scalar, minimal waiting time afte file has been modified in seconds (default = 10)
%
% See also FT_PREPROCESSING, LOAD

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

tic
memtic
st = dbstack;
functionname = st.name;
fprintf([functionname ' function started\n']);

% set defaults
cfg.timeout  = ft_getopt(cfg, 'timeout', Inf);
cfg.checkinterval  = ft_getopt(cfg, 'checkinterval', 1);
cfg.writerateapprox  = ft_getopt(cfg, 'writerateapprox', 1000000); % ~1 MByte
cfg.minmodwait  = ft_getopt(cfg, 'minmodwait', 10); % ~1 MByte



                                    
fprintf([functionname ' function initialized\n']);


while ~available([filename '.mat'],cfg.writerateapprox,cfg.minmodwait)
 pause(cfg.checkinterval)
 waited_seconds = toc;
 if waited_seconds > cfg.timeout
   ft_warning([functionname ' timeout after waiting for ' num2str(waited_seconds) 'seconds']);
   break
 end
end

disp(['... attempting to load ' filename '.mat' ' ...'])


fprintf([functionname ' function finished\n']);
waited_seconds = toc;
memtoc
end

function avail = available(filename,standard_write_rate_Byte_per_second,minmodwaitseconds)
present = (exist(filename)==2);
locked = true;
f = dir(filename);
if ~isempty(f) && present
    seconds_past_since_modification = (now()-f.datenum)*100000;
    seconds_to_write_approx = f.bytes/standard_write_rate_Byte_per_second;
    locked = (seconds_to_write_approx > seconds_past_since_modification) || (seconds_past_since_modification < minmodwaitseconds);
end
avail = present && ~locked;
end
