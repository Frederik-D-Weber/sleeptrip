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
function [numberECG] = getECGChannelNumbers(datachannel)
% datachannel = data.label
ecg_channel = ft_channelselection('*ECG*',datachannel);
numberECG = -1;
if ~isempty(ecg_channel)
    numberECG = find(strcmp(ecg_channel(1),datachannel));
else
    ecg_channel = ft_channelselection('*EKG*',datachannel);
    if ~isempty(ecg_channel)
        numberECG = find(strcmp(ecg_channel(1),datachannel));
    end
end

end