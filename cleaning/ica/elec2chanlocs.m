function chanlocs=elec2chanlocs(elec)

%convert Cartesian coordinates (different order/direction)
chanlocs=struct('labels',elec.label,'X',elec.chanpos(:,2),'Y',-1*elec.chanpos(:,1),'Z',elec.chanpos(:,3))';

%add all the other fields from Cartesian coors using eeglab func
chanlocs = convertlocs(chanlocs, 'cart2all');