%for pairs of start/end samples, merge any overlap
function intervalsOut=mergeOverlappingIntervals(intervalsIn)
nrow=size(intervalsIn,1);

% ----- find union of intervals
intervalsIn = sort(intervalsIn,2);          % make the pairs always increasing pairs
[intervalsIn, ind] = sort(intervalsIn(:));  % sorts all the input values, and keep track of how they moved around
c = [ones(1,nrow) -ones(1,nrow)]';          % this is a matrix that keeps track of when intervals start (1) and stop (-1)
c = c(ind);                                 % put in order of occurrence
csc = cumsum(c);                            %sum up starts (1) and stops (-1) , will be =0 at upper end of new interval(s)
irit = find(csc==0);                        % find index locations of 0 (ends of intervals)
ilef = [1; irit+1];                         % start of intervals index is at the very start (1) and one after all the other ends
ilef(end) = [];                             % no new interval starting at the very end

% intervalsOut matrix is start and end points of the new intervals
intervalsOut = [intervalsIn(ilef) intervalsIn(irit)];
end
