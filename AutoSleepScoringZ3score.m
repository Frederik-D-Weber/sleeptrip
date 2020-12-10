% Z3score V1 and V2 connector class
% Requires cfslib MATLAB library: https://github.com/amiyapatanaik/cfslib-MATLAB
% cfg.version     =  either 'z3scoreV1' or 'z3scoreV2'
% cfg.channel     =  a Nx1 number vector ordered in the with N = 4/5 channels in the data for 
%                    C3:A2
%                    C4:A1
%                    EOGl:A1/EOGl
%                    EOGr:A2/EOGr
%                    EMG (for V2 only)
%
% cfg.licensefile =  string, file path with a license file, haveing two
%                    lines. First line is username and second the api key,
%                    e.g,:
%                    example@email.com
%                    1234567890abcdefghijklmnopqrstuvwxyz12345678
%

% Copyright (C) 2017-, Frederik D. Weber
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
classdef AutoSleepScoringZ3score
    properties
        cfg
        scoring
        artifacts
        data
        isOperable
        isScored
    end
    methods
        function obj = AutoSleepScoringZ3score(cfg,data)
            obj.data = data;
%             if isempty(cfg)
                obj.cfg = [];
                obj.cfg.version = 'z3scoreV1';
                if isfield(cfg,'channel')
                    obj.cfg.channel = cfg.channel;
                end
                obj.cfg.autoSleepsCoreAlgoserverURL = 'https://z3score.com/api/v1';
                if isfield(cfg,'version')
                    if strcmp(cfg.version,'z3scoreV2')
                        obj.cfg.version = 'z3scoreV2';
                        obj.cfg.autoSleepsCoreAlgoserverURL = 'https://z3score.com/api/v2';
                    end
                end
                obj.cfg.email = 'example@email.com';
                obj.cfg.key = '1234567890abcdefghijklmnopqrstuvwxyz12345678';
                if isfield(cfg, 'licensefile')
                    if exist(cfg.licensefile) == 2
                        try
                            licfiletable = readtable(cfg.licensefile,'FileType','text','ReadVariableNames',false);
                            %obj.cfg.autoSleepsCoreAlgoserverURL = licfile{1};
                            obj.cfg.email = strtrim(licfiletable.Var1{1});
                            obj.cfg.key = strtrim(licfiletable.Var1{2});
                        catch err
                        	ft_error('Could not read the the Z3score license file as stated as %s check for its formating', cfg.license)
                        end
                    else
                        ft_error('Could not find the Z3score license file as stated as %s', cfg.license)
                    end
                else
%                     tempdirsplit = strsplit(pwd,filesep);
%                     for iPath = numel(tempdirsplit):1
%                         if exist([strjoin(tempdirsplit(1:iPath),filesep) filesep 'Z3Score_license.txt']) == 2
%                             try
%                                 licfile = read_mixed_csv([cfg.licdir filesep 'Z3Score_license.txt'],',');
%                                 %obj.cfg.autoSleepsCoreAlgoserverURL = licfile{1};
%                                 obj.cfg.email = licfile{1};
%                                 obj.cfg.key = licfile{2};
%                             catch err
%                             end
%                             break
%                         end
%                     end
                end
%             else
%                 obj.cfg = cfg;
%             end
            obj.cfg.epoch_duration_seconds = 30;
            
            obj.scoring = [];
            %obj.scoring.confidence = [];
            obj.artifacts = [];
            
            obj.isOperable = false;
            obj.isScored = false;
            
            [obj, ok] = obj.update_cfg();
            if ok
                obj = obj.update_channel();
                obj = obj.check();
            end
            
        end
        function [obj, ok] = update_cfg(obj)
            ok = true;
            prompt = {'Auto Sleep Score Algo Server URL', 'user email', 'API key'};
            title = 'Credentials Z3Score';
            dims = [1 35];
            definput = {obj.cfg.autoSleepsCoreAlgoserverURL, obj.cfg.email, obj.cfg.key};
            cfg_answers = inputdlg(prompt,title,dims,definput);
            if ~isempty(cfg_answers)
                obj.cfg.autoSleepsCoreAlgoserverURL = strtrim(cfg_answers{1});
                obj.cfg.email = strtrim(cfg_answers{2});
                obj.cfg.key = strtrim(cfg_answers{3});
            else
                disp(['Z3Score: abort did not update config!'])
                ok = false;
            end
        end
        function obj = update_channel(obj)
            if ~isfield(obj.cfg,'channel')
                obj.cfg.channel = [];
            else
                channelidx = obj.cfg.channel;
                obj.cfg.channel = [];
                obj.cfg.channel.C3_A2_idx = channelidx(1);
                obj.cfg.channel.C4_A1_idx = channelidx(2);
                obj.cfg.channel.EOGl_A1_idx = channelidx(3);
                obj.cfg.channel.EOGr_A2_idx = channelidx(4);
                if strcmp(obj.cfg.version, 'z3scoreV2')
                    obj.cfg.channel.EMG_idx = channelidx(5);
                end
            end
            %% choose channels
            if ~isfield(obj.cfg.channel,'C3_A2_idx')
                obj.cfg.channel.C3_A2_idx = 1;
            end
            if ~isfield(obj.cfg.channel,'C4_A1_idx')
                obj.cfg.channel.C4_A1_idx = 1;
            end
            if ~isfield(obj.cfg.channel,'EOGl_A1_idx')
                obj.cfg.channel.EOGl_A1_idx = 1;
            end
            if ~isfield(obj.cfg.channel,'EOGr_A2_idx')
                obj.cfg.channel.EOGr_A2_idx = 1;
            end
            switch obj.cfg.version
                case 'z3scoreV1'
                    channels_idx = [obj.cfg.channel.C3_A2_idx obj.cfg.channel.C4_A1_idx...
                        obj.cfg.channel.EOGl_A1_idx obj.cfg.channel.EOGr_A2_idx];
                    channel_list_request = {'C3:A2', 'C4:A1', 'EOGl:A1', 'EOGr:A2'};
                case 'z3scoreV2'
                    if ~isfield(obj.cfg.channel,'EMG_idx')
                        obj.cfg.channel.EMG_idx = 1;
                    end
                    channels_idx = [obj.cfg.channel.C3_A2_idx obj.cfg.channel.C4_A1_idx...
                        obj.cfg.channel.EOGl_A1_idx obj.cfg.channel.EOGr_A2_idx, obj.cfg.channel.EMG_idx];
                    channel_list_request = {'C3:A2', 'C4:A1', 'EOGl:A1', 'EOGr:A2', 'EMG'};
            end

            list_channel = obj.data.label;
            for iReqCh = 1:numel(channel_list_request)
                [selection,ok] = listdlg('PromptString',['Choose ' channel_list_request{iReqCh} ' channel'],...
                    'SelectionMode','single',...
                    'InitialValue',channels_idx(iReqCh),...
                    'ListString',list_channel);
                if ~ok
                    disp([channel_list_request{iReqCh} ' channel has not been updated'])
                    %return
                end
                switch iReqCh
                    case 1
                        obj.cfg.channel.C3_A2_idx = selection;
                    case 2
                        obj.cfg.channel.C4_A1_idx = selection;
                    case 3
                        obj.cfg.channel.EOGl_A1_idx = selection;
                    case 4
                        obj.cfg.channel.EOGr_A2_idx = selection;
                    case 5
                        obj.cfg.channel.EMG_idx = selection;
                    otherwise
                end
            end

        end
        function obj = update_data(obj,data)
            obj.data = data;
            obj = obj.update_channel();
        end
        function message = status(obj)
            %check license
            try
                %response = loadjson(urlread([obj.cfg.autoSleepsCoreAlgoserverURL '/check' '?' 'email=' obj.cfg.email '&' 'key=' obj.cfg.key]));

                response = loadjson(urlreadpost([obj.cfg.autoSleepsCoreAlgoserverURL '/check'],...
                    {'email', obj.cfg.email, 'key', obj.cfg.key}));
            catch
                message = 'Z3Score: scoring server is unreachable';
                disp(message);
                return
            end
            
            if response.status == 0,
                message = 'Z3Score: License check failed. Check username and API-key';
                disp(message);
                disp(['Z3Score: Error message: ' response.message])
                return
            end
            
            %disp(response.message);
            
            num_epochs_requested = floor(size(obj.data.trial{1},2)/obj.data.fsample/obj.cfg.epoch_duration_seconds);

            message = [sprintf('Z3Score: Server reachable\nAPI Credits limit left (hourly): %d\n Epoch limit left (epochs): %d\n ... next request %d epochs', response.call_limit, response.epoch_limit, num_epochs_requested)];
            
        end
        function obj = check(obj)
            obj.isOperable = false;
            %check license
            try
                response = loadjson(urlreadpost([obj.cfg.autoSleepsCoreAlgoserverURL '/check'],...
                    {'email',obj.cfg.email,'key',obj.cfg.key}));
            catch
                disp('Z3Score: scoring server is unreachable');
                return
            end
            
            if response.status == 0,
                disp('Z3Score: License check failed');
                disp(['Z3Score: Error message: ' response.message])
                return
            end
            
            %disp(response.message);
            
            num_epochs_requested = floor(size(obj.data.trial{1},2)/obj.data.fsample/obj.cfg.epoch_duration_seconds);

            if ~(response.call_limit >= 1) && ~((response.epoch_limit-num_epochs_requested) > 0)
                msgbox('Limit reached' ,sprintf('Z3Score: API Call limit left (hourly): %d\n Epoch limit left (daily): %d\n ... but requested %d', response.call_limit, response.epoch_limit, num_epochs_requested), 'modal');
                return
            end
            
            
            if isempty(obj.cfg.channel.C3_A2_idx) &&...
                isempty(obj.cfg.channel.C4_A1_idx) &&...
                isempty(obj.cfg.channel.EOGl_A1_idx) &&...
                isempty(obj.cfg.channel.EOGr_A2_idx)
                if strcmp(obj.cfg.version,'z3scoreV2')
                    if isempty(obj.cfg.channel.EMG_idx)
                        disp('Z3Score: data does not provide the 5 needed channels')
                    end
                else
                %size(obj.data.trial{1},1) ~= 4
                disp('Z3Score: data does not provide the 4 needed channels')
                end
                return
            end
            obj.isOperable = true;
        end
        function obj = score(obj)
            if obj.isOperable
                if ~isfield(obj.cfg,'stream')
                    obj.cfg.stream = [];
                end
                if isempty(obj.cfg.stream)
                    %Find out sampling rate
                    samplingRate = obj.data.fsample;
                    %Construct raw data from selected channels
                    switch obj.cfg.version
                        case 'z3scoreV1'
                            channels_idx = [obj.cfg.channel.C3_A2_idx obj.cfg.channel.C4_A1_idx...
                                obj.cfg.channel.EOGl_A1_idx obj.cfg.channel.EOGr_A2_idx];
                        case 'z3scoreV2'
                            channels_idx = [obj.cfg.channel.C3_A2_idx obj.cfg.channel.C4_A1_idx...
                                obj.cfg.channel.EOGl_A1_idx obj.cfg.channel.EOGr_A2_idx obj.cfg.channel.EMG_idx];
                    end
  
                    EEGData = obj.data.trial{1}(channels_idx,:);
                    
                    num_epochs = floor(size(EEGData,2)/samplingRate/obj.cfg.epoch_duration_seconds);
                    
                    
                    
                    %Convert raw stream to a CFS and write to a file
                    disp('Z3Score: put data to CFS stream ......... this may take some time depending on your connection!');
                    %    tic;
                    %Convert raw PSG stream to CFS stream
                    switch obj.cfg.version
                        case 'z3scoreV1'
                           [stream, status, message, quality] = streamCFS(EEGData, samplingRate);
                           disp(['Z3Score data message: ' message]);
                           obj.cfg.stream = stream;
                        case 'z3scoreV2'
                           [stream, status, message, quality] = streamCFS_V2(EEGData(1,:)', EEGData(2,:)', EEGData(3,:)', EEGData(4,:)', EEGData(5,:)', samplingRate, samplingRate, samplingRate);
                           disp(['Z3Score data message: ' message]);
                           obj.cfg.stream = stream;
                    end

                    %obj.cfg.stream = streamCFS_V2(signalCell{C3N}, signalCell{C4N}, signalCell{ELN}, signalCell{ERN}, signalCell{EM}, samplingRateEEG, samplingRateEOG, samplingRateEMG);
                    %obj.cfg.stream = streamCFS_V2(EEGData(1,:), EEGData(2,:), EEGData(3,:), EEGData(4,:), EEGData(5,:), samplingRate, samplingRate, samplingRate);

                    EEGData = [];
                    %fileID = fopen('test.cfs','w');
                    %fwrite(fileID,stream,'*uint8','ieee-le');
                    %fclose(fileID);
                    %t = toc;
                    %fprintf('Time taken %.3f seconds\n',t);
                end
                
                disp('Z3Score: fetch scoring info from server');
     
                try
                    response = loadjson(urlreadpost([obj.cfg.autoSleepsCoreAlgoserverURL  '/score'], ...
                        {'email',obj.cfg.email,'key',obj.cfg.key,'file',obj.cfg.stream}));
                catch
                    disp('Z3Score: Scoring server is unreachable');
                    return
                end
                
                if response.status == 0,
                    disp('Z3Score: Error scoring data.');
                    disp(['Z3Score: Error message:' response.message])
                    return
                end
                
                %fprintf('Time taken %.3f seconds.\nAPI calls left (hourly limits): %d, Epochs left (daily limits): %d \n',t, response.calls_left, response.epochs_left);
                %Automatic sleep scores
                scores = response.message;
                obj.scoring = [];
                obj.scoring.epochlength = 30;
                obj.scoring.epochlength_2 = 5;
                obj.scoring.epochs = scores(:,1)';
                obj.scoring.epochs_2 = scores(:,1)';
                obj.scoring.excluded = logical(zeros(1,numel(scores(:,1))));
                obj.scoring.confidence = scores(:,2)';
                obj.scoring.confidence = linear_scaling_normalization(obj.scoring.confidence, 0, 10, 0, 1);

                %artifacts are in 5-seconds scoring, only 7 are
                %unscorable/artifact epochs
                artifacts_5sec_blocks = response.artifact;
                obj.scoring.epochs_2 = artifacts_5sec_blocks';
                artifacts_5sec_blocks_index = (artifacts_5sec_blocks == 7);
                
                epoch_length_seconds = 5; % in seconds
                temp_starting_seconds_in_data = (1:floor(obj.data.sampleinfo(2)/(obj.data.fsample*epoch_length_seconds)))-1;
                epoch_begsample = obj.data.fsample*epoch_length_seconds*temp_starting_seconds_in_data+1;
                epoch_endsample = obj.data.fsample*epoch_length_seconds*(temp_starting_seconds_in_data+1);
                obj.artifacts = [epoch_begsample' epoch_endsample'];
                obj.artifacts = obj.artifacts(artifacts_5sec_blocks_index,:);
                
                %artifacts MA epochs mark
                for iEpoch = 1:numel(obj.scoring.epochs)
                    obj.scoring.excluded(iEpoch) = any(artifacts_5sec_blocks_index(((iEpoch-1)*6+1):(iEpoch*6)))*3;
                end
                
                obj.scoring.ori = [];
                obj.scoring.ori.epochs = obj.scoring.epochs;
                obj.scoring.ori.excluded = obj.scoring.excluded;
                obj.scoring.ori.scores = scores;
                obj.scoring.ori.response = response;
                obj.scoring.ori.artifacts = obj.artifacts;
                
                obj.scoring.epochs = arrayfun(@num2str,obj.scoring.epochs,'UniformOutput',false);
                
                scoremap = [];
                scoremap.labelold  = {'0', '1',  '2',  '3',  '4',  '5'};
                scoremap.labelnew  = {'W', 'N1', 'N2', 'N3', 'N3', 'R'};
                scoremap.unknown   = '?';
                
                obj.scoring.label = unique(scoremap.labelnew)';
                obj.scoring.label = ({obj.scoring.label{:}})';%assure the vertical orientation
                
                obj.scoring.standard = 'custom';
                
                cfg_sc = [];
                cfg_sc.scoremap = scoremap;
                cfg_sc.to = 'aasm';
                obj.scoring = st_scoringconvert(cfg_sc, obj.scoring);
                
                obj.scoring.ori.scoremap  = scoremap;
                obj.scoring.ori.scoringfile   = ['Z3score-API' ' ' obj.cfg.version];
                obj.scoring.ori.scoringformat   = 'aasm';
                
                obj.scoring.ori.cfg = obj.cfg;
                
                obj.scoring.ori.cfg = rmfield(obj.scoring.ori.cfg,'key'); %delete the private license key info
                disp('Z3Score: Done');
            else
                obj.check(obj);
            end
            function newvalue = linear_scaling_normalization(toNorm, minInterval, maxInterval, minNormInterval, maxNormInterval)
                newvalue =  ((toNorm - minInterval) / (maxInterval - minInterval)) * (maxNormInterval - minNormInterval) + minNormInterval;
            end
        end
    end
end




