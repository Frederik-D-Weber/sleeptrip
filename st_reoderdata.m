function [data] = st_reoderdata(cfg, data)

% ST_REORDERDATA applies filter strategy to the data.
%
% Use as
%   [data] = st_reoderdata(cfg, data)
%
% Available configuration parameters are:
%   cfg.order  = either number vector Nx1 or string like
%                'alphabetical'
%                'reverse_alphabetical'
%                'reverse'
%                'random'
%
% See also FT_PREPROCESSING, FT_APPLY_MONTAGE

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
cfg.order  = ft_getopt(cfg, 'order', 1:numel(data.label));

fprintf([functionname ' function initialized\n']);


if isnumeric(cfg.order)
    if size(cfg.order,1) == numel(data.label)
        curr_chanOrder = cfg.order;
    else
        ft_error('cfg.order of length %d must be Nx1 number vector of with N = %d.',numel(cfg.order), numel(data.label))
    end
else
    switch cfg.order
        case 'reverse'
            curr_chanOrder = numel(data.label):-1:1;
        case 'alphabetical'
            [tempchlab curr_chanOrder] = sort(data.label);
        case 'reverse_alphabetical'
            [tempchlab curr_chanOrder] = sort(data.label);
            curr_chanOrder = flip(curr_chanOrder);
        case 'random'
            curr_chanOrder = randperm(numel(data.label));
        case 'scoring'
            curr_chanOrder = randperm(numel(data.label));
    end
end

curr_chanOrder = curr_chanOrder(:);

% data_reord = {};
% iChannel = 1;
% for iChannelbyOrder = curr_chanOrder'
%     cfg_tmp = [];
%     cfg_tmp.channel = iChannelbyOrder;
%     data_reord{iChannel} = ft_selectdata(cfg_tmp,data);
%     iChannel = iChannel + 1;
% end
% 
% if length(data_reord) > 1
%     data = ft_appenddata([],data_reord{:});
% else
%     data = data_reord{:};
% end

nTr = numel(data.trial);
for iTr  = 1:nTr
    data.trial{iTr} = data.trial{iTr}(curr_chanOrder,:);
end

data.label = data.label(curr_chanOrder);

fprintf([functionname ' function finished\n']);
toc(ttic)
memtoc(mtic)
end