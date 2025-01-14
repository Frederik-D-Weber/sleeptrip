function [res_sp_train_channel, res_sp_train_event] = st_spindle_train(cfg, res_event_sp)

% ST_SPINDLE_TRAIN finds spindle trains from the results of st_spindles.
%
% Use as:
%   [res_sp_train_channel, res_sp_train_event] = st_spindle_train(cfg, res_event_sp)
%
% Optional configuration parameters are:
%
% cfg.time_inter_SP_cutoff = Maximal number of seconds between two consecutive spindles on the same EEG channel to be considered as grouped (default = 6)
% cfg.min_nbr_SP_in_train  = Minimal number of grouped spindles to define a train (default = 2)
%
% The DEFAULT parameters are based on Boutin and Doyon (2020), defining
% spindle train as two or more consecutive and electrode-specific spindle 
% events interspaced by less than or equal to 6 s. This definition is the
% most widely used definition (see for instance Solano et al., 2022;
% Champetier et al., 2023 or Boutin et al., 2024)

% Copyright (C) 2024-, Pierre Champeter
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

%% Display function name when user runs the function

functionname = getfunctionname();
fprintf([functionname ' function started\n']);


%% Constant parameters

% 1) time_inter_SP_cutoff
if isfield(cfg, 'time_inter_SP_cutoff')
    time_inter_SP_cutoff = cfg.time_inter_SP_cutoff;
else
    time_inter_SP_cutoff = 6;
end

% 2) min_nbr_SP_in_train
if isfield(cfg, 'min_nbr_SP_in_train')
    min_nbr_SP_in_train = cfg.min_nbr_SP_in_train;
else
    min_nbr_SP_in_train = 2;
end


%% Define variables

% 1) Create main outputs variables
    % res_sp_train_channel
res_sp_train_channel = [];
res_sp_train_channel.ori = 'st_spindle_train';
res_sp_train_channel.type = 'spindles_train_channel';
res_sp_train_channel.cfg = cfg;
res_sp_train_channel.table = table();

    % res_sp_train_event
res_sp_train_event = [];
res_sp_train_event.ori = 'st_spindle_train';
res_sp_train_event.type = 'spindles_train_event';
res_sp_train_event.cfg = cfg;
res_sp_train_event.table = table();

% 2)Identify EEG channel names
channel_list = unique(res_event_sp.table.channel);

% 3) Initialize variable before entering into for loop
time_sec_bw_SP = table();


%% For loop start (1 iteration per EEG channel)

for nE = 1:size(channel_list,1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        I - res_sp_train_channel (average values for each channel)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% I- 1) Select spindle events of the electrode iter (nE)

    row_nE_iter = ismember(res_event_sp.table.channel, channel_list(nE));

    sp_nE_iter =  res_event_sp.table(row_nE_iter,:);


    %% I- 2) Calculate the time interval between each spindle event

    time_sec_bw_SP_nE_iter = [];
    for sp_iter = 1:size(sp_nE_iter,1)-1
        time_sec_bw_SP_nE_iter = [time_sec_bw_SP_nE_iter,  sp_nE_iter.seconds_begin(sp_iter+1) - sp_nE_iter.seconds_end(sp_iter)];
    end


    %% I- 3) Determine fast SP trains characteristics

    % 1) Calculate the length of all the potential SP groups (can be equal to 1)
    length_sp_candidate_train = [];
    length_iter = 1;
    for i = 1:size(time_sec_bw_SP_nE_iter, 2)
        if time_sec_bw_SP_nE_iter(i) <= time_inter_SP_cutoff && i~= size(time_sec_bw_SP_nE_iter, 2)
            length_iter = length_iter +1;
        elseif time_sec_bw_SP_nE_iter(i) > time_inter_SP_cutoff
            length_sp_candidate_train = [length_sp_candidate_train, length_iter];
            length_iter = 1;
        elseif time_sec_bw_SP_nE_iter(i) <= time_inter_SP_cutoff && i== size(time_sec_bw_SP_nE_iter, 2)
            length_iter = length_iter +1;
            length_sp_candidate_train = [length_sp_candidate_train, length_iter];
        end
    end

    % 2) Identify spindle trains (= the candidate trains that contain at least min_nbr_SP_in_train spindles)
    index_sp_train = find(length_sp_candidate_train >= min_nbr_SP_in_train);

    % 3) Determine the fast SP trains characteristics
    table_iter = table();


    sleep_stages = res_event_sp.table.used_stages_for_detection(1);
    channel = channel_list(nE);
    nb_sp_trains = size(index_sp_train,2);
    mean_nb_sp_per_sp_train = mean(length_sp_candidate_train(index_sp_train));
    median_nb_sp_per_sp_train = median(length_sp_candidate_train(index_sp_train));
    max_nb_sp_per_sp_train = max(length_sp_candidate_train(index_sp_train));
    min_nb_sp_per_sp_train = min(length_sp_candidate_train(index_sp_train));
    nb_sp_in_trains = sum(length_sp_candidate_train(find(length_sp_candidate_train >= min_nbr_SP_in_train)));
    prop_sp_in_trains = 100 * nb_sp_in_trains / sum(length_sp_candidate_train);


    %% I- 4) Gather results in res_sp_train_channel

    table_iter = table(channel, sleep_stages, nb_sp_trains, ...
        mean_nb_sp_per_sp_train, median_nb_sp_per_sp_train, max_nb_sp_per_sp_train, min_nb_sp_per_sp_train, ...
        nb_sp_in_trains, prop_sp_in_trains);

    res_sp_train_channel.table = [res_sp_train_channel.table; table_iter];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       II - res_sp_train_event
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% II - 1) Find the begin and end (i.e., 1st spindle and last spindle) of each train candidate

    index_sp_train_candidate_begin_end = [];

    begin_iteratif = 1;
    end_iteratif = 1;

    for i = 1:size(time_sec_bw_SP_nE_iter, 2)
        if time_sec_bw_SP_nE_iter(i) <= time_inter_SP_cutoff && i~= size(time_sec_bw_SP_nE_iter, 2)
            end_iteratif = end_iteratif + 1;
        elseif time_sec_bw_SP_nE_iter(i) > time_inter_SP_cutoff && i~= size(time_sec_bw_SP_nE_iter, 2)
            index_sp_train_candidate_begin_end = [index_sp_train_candidate_begin_end; begin_iteratif, end_iteratif];
            begin_iteratif = i + 1;
            end_iteratif = i + 1;
            
            %Add the last SP train if it is the end of the array
        elseif time_sec_bw_SP_nE_iter(i) > time_inter_SP_cutoff && i== size(time_sec_bw_SP_nE_iter, 2)
            index_sp_train_candidate_begin_end = [index_sp_train_candidate_begin_end; begin_iteratif, end_iteratif];
            index_sp_train_candidate_begin_end = [index_sp_train_candidate_begin_end; i+1, i+1];
            
            %Add the last isolated SP if it is the end of the array
        elseif time_sec_bw_SP_nE_iter(i) <= time_inter_SP_cutoff && i== size(time_sec_bw_SP_nE_iter, 2)
            index_sp_train_candidate_begin_end = [index_sp_train_candidate_begin_end; begin_iteratif, i+1];
        end
    end


%% II - 2) Compute characteristics for each candidate spindle train

    % 1) Initialize variables 
    Nb_sp_in_the_train_candidate             = [];
    Mean_sp_duration_in_the_train_candidate  = [];
    Mean_sp_ampl_in_the_train_candidate      = [];
    Mean_sp_freq_in_the_train_candidate      = [];
    Begin_train_candidate_sec                = [];
    End_train_candidate_sec                  = [];
    Begin_sec_each_sp_in_the_train_candidate = [];
    End_sec_each_sp_in_the_train_candidate   = [];
    Train_candidate_duration_sec             = []; 

    % 2) For all spindle train candidate
    for row = 1:size(index_sp_train_candidate_begin_end,1)
        number_fast_SP_in_train_iter = index_sp_train_candidate_begin_end(row,2) - index_sp_train_candidate_begin_end(row,1) + 1;
        mean_fast_SP_duration_of_each_train_iter = mean(sp_nE_iter.duration_seconds(index_sp_train_candidate_begin_end(row,1):index_sp_train_candidate_begin_end(row,2)) );
        mean_fast_SP_ampl_of_each_train_iter = mean(sp_nE_iter.amplitude_peak2trough_max(index_sp_train_candidate_begin_end(row,1):index_sp_train_candidate_begin_end(row,2)) );
        mean_fast_SP_freq_of_each_train_iter = mean(sp_nE_iter.frequency_by_mean_pk_trgh_cnt_per_dur(index_sp_train_candidate_begin_end(row,1):index_sp_train_candidate_begin_end(row,2)) );

        Nb_sp_in_the_train_candidate = [Nb_sp_in_the_train_candidate; number_fast_SP_in_train_iter];
        Mean_sp_duration_in_the_train_candidate = [Mean_sp_duration_in_the_train_candidate; mean_fast_SP_duration_of_each_train_iter];
        Mean_sp_ampl_in_the_train_candidate = [Mean_sp_ampl_in_the_train_candidate; mean_fast_SP_ampl_of_each_train_iter];
        Mean_sp_freq_in_the_train_candidate = [Mean_sp_freq_in_the_train_candidate; mean_fast_SP_freq_of_each_train_iter];

        % Add begin and end of the train candidate
            % Begin train = begin of the 1st SP of the train
        begin_sec_iter = sp_nE_iter.seconds_begin(index_sp_train_candidate_begin_end(row,1));
            % End train = end of the last SP of the train
        end_sec_iter = sp_nE_iter.seconds_end(index_sp_train_candidate_begin_end(row,2));

        Begin_train_candidate_sec = [Begin_train_candidate_sec; begin_sec_iter];
        End_train_candidate_sec = [End_train_candidate_sec; end_sec_iter];
        Middle_train_candidate_sec = (Begin_train_candidate_sec + End_train_candidate_sec)/2;
   
        if number_fast_SP_in_train_iter == 1
            Train_candidate_duration_sec = [Train_candidate_duration_sec; mean_fast_SP_duration_of_each_train_iter];
        else
            Train_candidate_duration_sec = [Train_candidate_duration_sec; End_train_candidate_sec(row) - Begin_train_candidate_sec(row)];
        end
        % Add begin and end of each SP in each train
        Begin_sec_each_sp_in_the_train_candidate = [Begin_sec_each_sp_in_the_train_candidate; {sp_nE_iter.seconds_begin(index_sp_train_candidate_begin_end(row,1):index_sp_train_candidate_begin_end(row,2))}];
        End_sec_each_sp_in_the_train_candidate = [End_sec_each_sp_in_the_train_candidate ;{sp_nE_iter.seconds_end(index_sp_train_candidate_begin_end(row,1):index_sp_train_candidate_begin_end(row,2))}];
    end

    % 3) Put the data in the table
    channel = sp_nE_iter.channel(index_sp_train_candidate_begin_end(:,1));
    sleep_stages = sp_nE_iter.used_stages_for_detection(index_sp_train_candidate_begin_end(:,1));

    iter_electr_trains = table(...
        channel,...
        sleep_stages,...
        Nb_sp_in_the_train_candidate, ...
        Mean_sp_duration_in_the_train_candidate, ...
        Mean_sp_ampl_in_the_train_candidate, ...
        Mean_sp_freq_in_the_train_candidate, ...
        Begin_train_candidate_sec, ...
        End_train_candidate_sec, ...
        Middle_train_candidate_sec, ...
        Train_candidate_duration_sec, ...
        Begin_sec_each_sp_in_the_train_candidate, ...
        End_sec_each_sp_in_the_train_candidate ...
        );

    % 4) Gather the results of all the electrodes in a single table
    res_sp_train_event.table = [res_sp_train_event.table; iter_electr_trains];

end


end