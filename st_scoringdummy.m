function scoring=st_scoringdummy(cfg)

%cfg.epochlength=ft_getopt(cfg, 'epochlength', 30);
cfg.label=ft_getopt(cfg, 'label', {'?'});
cfg.epochnumber=ft_getopt(cfg, 'epochnumber', 1000);

%simple stage table
tableScoring=table(repmat(cfg.label,[cfg.epochnumber,1]),'VariableNames',{'Stage'});

%set up scoremap
scoremap = [];
scoremap.labelold  = {'W', 'N1', 'N2', 'N3', 'R'};
scoremap.labelnew  = {'W', 'N1', 'N2', 'N3', 'R'};
scoremap.unknown   = '?';

%read the scoring
cfg = [];
cfg.standard='custom';
cfg.to='aasm';
cfg.scoremap=scoremap;

scoring=st_read_scoring(cfg,tableScoring);