function [begins, ends] = getBeginsAndCorrespondingEndsThresholdCrossing(data,threshold,isRisingAboveNotFallingBelow)

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

    tx = getThresholdCrossing(data,threshold,isRisingAboveNotFallingBelow);
    if length(tx) >= 2
        if mod(length(tx),2) ~= 0
            tx = tx(1:end-1);
        end
        begins = tx(1:2:end);
        ends = tx(2:2:end);
    else
        begins = [];
        ends = [];
    end
end


function tx = getThresholdCrossing(data,threshold,isRisingAboveNotFallingBelow)
% get indices of threshold crossing in data starting with a rising/falling edge,
%i.e. from below to above threshold
    data = data - threshold;
    temp = data(1:end-1) .* data(2:end); % scalar multiplication of two arrays. Second array is shifted to left by 1.
    p = find(temp<0);
    if length(p) >= 2
        if isRisingAboveNotFallingBelow
            if(data(p(1))<0) % transition from below to above.
                tx = p(1:1:end);
            else             % transition from above to below.
                tx = p(2:1:end);
            end
        else
            if(data(p(1))<0) % transition from above to below.
                tx = p(2:1:end);
            else             % transition from below to above.
                tx = p(1:1:end);
            end
        end
        tx=tx';
    else
        tx = [];
    end
    
end