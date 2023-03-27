
%for pairs of start/end samples, merge if gap between pairs smaller than mergeWindow
function intervalsOut=mergeCloseIntervals(intervalsIn,mergeWindow)
intervalsIn = sort(intervalsIn,2);          % make the pairs increasing (each pair on row)
intervalsIn = sort(intervalsIn,1);          %also ensure we start with earliest (possibly negative) pair

%sort based on interval start
[~, index] = sort(intervalsIn(:,1),1);
intervalsIn   = intervalsIn(index,:);

%
mergeStatusChange=diff([0;[intervalsIn(2:end,1)-intervalsIn(1:end-1,2);nan]<=mergeWindow]);
mergeInds=[find(mergeStatusChange==1) find(mergeStatusChange==-1)];

if ~isempty(mergeInds)
    newIntNum=size(mergeInds,1);
    newIntervals=zeros(newIntNum,2);
    for row_i=1:newIntNum
        newIntervals(row_i,:)=[intervalsIn(mergeInds(row_i,1),1) intervalsIn(mergeInds(row_i,2),2)];
    end

    mergeIndsCell=num2cell(mergeInds);
    removeInds=cellfun(@colon,mergeIndsCell(:,1),mergeIndsCell(:,2),'uni',false);
    intervalsIn(unique(cat(2,removeInds{:})),:)=[];

    intervalsOut=[intervalsIn;newIntervals];

    %sort again
    [~, index] = sort(intervalsOut( :,1), 1);
    intervalsOut   = intervalsOut(index,:);
else
    intervalsOut=intervalsIn;
end