function [memoryavailGB, avail, waited_seconds] = st_wait_memory(cfg, memoryneededGB)

% ST_WAIT_EXISTS pauses program until filename exists
%
% Use as
%   [memoryneededGB]                        = st_wait_exists(cfg, memoryneededGB)
%   [memoryneededGB, avail]                 = st_wait_exists(cfg, memoryneededGB)
%   [memoryneededGB, avail, waited_seconds] = st_wait_exists(cfg, memoryneededGB)
%
%   cfg can be empty (i.e., cfg = [];)
%
% Optional configuration parameters are:
%
%   cfg.timeout         = scalar, time in seconds to wait at most (default = Inf)
%   cfg.checkinterval   = scalar, interval in seconds to check (default = 1)
%   cfg.minmemoryneededGB = scalar, minimal GByte to remain free in addition as a buffer (default = 0.5)
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

ttic = tic;
mtic = memtic;
functionname = getfunctionname();
fprintf([functionname ' function started\n']);

% set defaults
cfg.timeout  = ft_getopt(cfg, 'timeout', Inf);
cfg.checkinterval  = ft_getopt(cfg, 'checkinterval', 1);
cfg.minmemoryneededGB  = ft_getopt(cfg, 'minmemoryneededGB', 0.5); % 0.5 GByte

                                    
fprintf([functionname ' function initialized\n']);

fprintf(['waiting for %f GByte to be available with minimum buffer of %f (%f in total)' '\n'],memoryneededGB,cfg.minmemoryneededGB,memoryneededGB+cfg.minmemoryneededGB);


actualneededmemory_inBytes = (memoryneededGB+cfg.minmemoryneededGB)*1024*1024*1024;
[userview,systemview] = memory;
present = userview.MemAvailableAllArrays;

avail = present >= actualneededmemory_inBytes;

while ~avail
 pause(cfg.checkinterval)
 waited_seconds = toc(ttic);
 if waited_seconds > cfg.timeout
   ft_warning([functionname ' timeout after waiting for ' num2str(waited_seconds) 'seconds']);
   break
 end
 
[userview, systemview] = memory;
present = userview.MemAvailableAllArrays;
avail = present >= actualneededmemory_inBytes;

end

memoryavailGB = userview.MemAvailableAllArrays/1024/1024/1024;

waited_seconds = toc(ttic);
fprintf(['waited for ' num2str(waited_seconds) ' seconds ' 'for ' memoryneededGB ' to be finally available'  '\n']);

fprintf([functionname ' function finished\n']);
memtoc(mtic)
end

 
