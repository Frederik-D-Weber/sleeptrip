function data=st_reref(cfg)

ttic = tic;
mtic = memtic;
functionname = getfunctionname();
fprintf([functionname ' function started\n']);

% these are used by the ft_preamble/ft_postamble function and scripts
ft_revision = '$Id$';
ft_nargin   = nargin;
ft_nargout  = nargout;

% do the general setup of the function
st_defaults
ft_preamble init
ft_preamble debug
ft_preamble loadvar data
ft_preamble provenance data
ft_preamble trackconfig

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
    return
end

%-----core function start---
%---input checks and defaults----
ft_checkconfig(cfg,'required',{'data','reflabel'});
cfg.removelabel  = ft_getopt(cfg, 'removelabel', {}); % default, nothing to remove

fprintf([functionname ' function initialized\n'])


elecLabels=cfg.data.label;
numChans=size(elecLabels,1);

%locate channels to become new reference
refLabels=cfg.reflabel;
refInds=[];
for refLabel_i=1:length(refLabels)
    refInd=find(strcmp(elecLabels,refLabels{refLabel_i}));

    if isempty(refInd)
        fprintf('reference label not found\n')
        rereferenceStruct=struct;
        return
    end

    refInds=[refInds refInd];
end

%locate channels to remove from rereferenced data
removeLabels=cfg.removelabel;
removeInds=[];
for removeLabel_i=1:length(removeLabels)
    removeInd=find(strcmp(elecLabels,removeLabels{removeLabel_i}));

    if isempty(removeInd)
        fprintf('remove label not found\n')
        rereferenceStruct=struct;
        return
    end

    removeInds=[removeInds removeInd];
end

%add ref chans to removal chans
removeInds=[refInds removeInds];

refMat=zeros(numChans,numChans);
refMat(1:numChans+1:end)=1; %ones for diagonal
refMat(:,refInds)=-1/length(refInds); %e.g., -0.5 for mastoids
refMat(removeInds,:)=[]; %remove to-be-esxcluded channels

%build Fieldtrip reref struct
rereferenceStruct.tra=refMat;
rereferenceStruct.labelold = elecLabels;
rereferenceStruct.labelnew = elecLabels; %keep the labels...
rereferenceStruct.labelnew(removeInds,:)=[]; %..%but remove to-be-removed chans

%rereference
data=ft_apply_montage(cfg.data,rereferenceStruct);

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
%ft_postamble trackconfig
% ft_postamble previous data
% data = datanew
% ft_postamble provenance data
% ft_postamble history    data
% ft_postamble savevar    data

fprintf([functionname ' function finished\n']);
toc(ttic)
memtoc(mtic)

