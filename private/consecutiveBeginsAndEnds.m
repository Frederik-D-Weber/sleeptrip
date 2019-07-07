% find begins and end samples in a Nx2 vector and concatenate the beginnin and ends to a shorter Nx2 if possible

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

function [begins, ends] =  consecutiveBeginsAndEnds(x,maxConsecutiveDiff)
    nx = length(x);
    iter = 1;
    begins(iter) = x(1);
    ends = [];
    if nx > 1
        for i = 2:nx
            if x(i) > (x(i-1) + maxConsecutiveDiff)
                ends(iter) = x(i-1);
                iter = iter + 1;
                begins(iter) = x(i);
            end
        end
    end
    if length(ends) < length(begins)
        ends(iter) = x(nx);
    end
end