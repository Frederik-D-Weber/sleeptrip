%smoothing,  with sliding window of wd data points, i.e. the center sample +-floor(wd/2), 
%therefore wd is recomended to be an odd positive integer and smaller than points in the data. 
%At the boundaries the mean is smoothed using a reduced window size that shrinks to the maximum samples left over to the left or right of the center sample.
%smooth(mean(tfa.powspctrm,1),round(SmoothFreqWindowSize/(1/SegmentLength)),'lowess');

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

function r = smoothwd(x,wd)
if (mod(wd,2) == 0) || (wd < 3) || (wd > length(x))
    error('the window size wd is supposed to be an positive odd number greater than 2, but smaller than the number of data points x!');
end

l = length(x);
r = zeros(1,l);
lastSample = l - wd + 1;
tempAddISmpl = wd - 1;
tempWindows = zeros(lastSample,wd);

for iSmpl = 1:lastSample
    tempWindows(iSmpl,:) = x(iSmpl:(iSmpl + tempAddISmpl));
end

firstSmplstimeWndw = (wd+1)/2 ;
lastSmplstimeWndw = l - ((wd-1)/2);
r(firstSmplstimeWndw:lastSmplstimeWndw) = mean(tempWindows,2);

r(1:(firstSmplstimeWndw-1)) = smoothLeftBorder(x(1:((firstSmplstimeWndw*2)-1)),wd);
r((lastSmplstimeWndw+1):end) = smoothRightBorder(x(((lastSmplstimeWndw-firstSmplstimeWndw)+1):end),wd);

end


function r = smoothRightBorder(x,wd)
bu = floor(wd/2);
r = zeros(1,bu);
for i=1:bu
    r(i) = mean(x((2*i-1):(2*bu-1)));
end
end


function r = smoothLeftBorder(x,wd)
    r = smoothRightBorder(x(end:-1:1),wd);
    r = r(end:-1:1);
end