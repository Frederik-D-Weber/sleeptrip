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
function [numberEEG numberEEG_frontal numberEEG_occipital numberEOG numberEMG numberECG] = getScoringChannelNumbers(datachannel)
% datachannel = data.label

numberEEG = -1;
numberEEG_frontal = -1;
numberEEG_occipital = -1;
numberEOG = -1;
numberEMG = -1;
numberECG = -1;
eeg_channel = ft_channelselection('*C4*',datachannel);
eeg_frontal_channel = ft_channelselection('*F4*',datachannel);
eeg_occipital_channel = ft_channelselection('*O2*',datachannel);
eog_channel = ft_channelselection('*EOG*',datachannel);
emg_channel = ft_channelselection('*EMG*',datachannel);
ecg_channel = ft_channelselection('*ECG*',datachannel);
ppg_channel = ft_channelselection({'OXY_IR_AC','*PPG*'},datachannel);


if ~isempty(eeg_channel)
    numberEEG = find(strcmp(eeg_channel(1),datachannel));
else
    if numel(datachannel) > 2
        numberEEG = 2;
    else
        numberEEG = 1;
    end
end

if ~isempty(eeg_frontal_channel)
    numberEEG_frontal = find(strcmp(eeg_frontal_channel(1),datachannel));
else
    if numel(datachannel) > 2
        numberEEG_frontal = 2;
    else
        numberEEG_frontal = 1;
    end
end

if ~isempty(eeg_occipital_channel)
    numberEEG_occipital = find(strcmp(eeg_occipital_channel(1),datachannel));
else
    if numel(datachannel) > 2
        numberEEG_occipital = 2;
    else
        numberEEG_occipital = 1;
    end
end

if ~isempty(eog_channel)
    numberEOG = find(strcmp(eog_channel(1),datachannel));
else
    numberEOG = 1;
end

if ~isempty(emg_channel)
    numberEMG = find(strcmp(emg_channel(1),datachannel));
else
    numberEMG = numel(datachannel);
end


if ~isempty(ecg_channel)
    numberECG = find(strcmp(ecg_channel(1),datachannel));
else
    ecg_channel = ft_channelselection('*EKG*',datachannel);
    if ~isempty(ecg_channel)
        numberECG = find(strcmp(ecg_channel(1),datachannel));
    else
        % take the PPG channel instead
        if ~isempty(ppg_channel)
        	numberECG = find(strcmp(ppg_channel(1),datachannel));
        end
    end
end

end

