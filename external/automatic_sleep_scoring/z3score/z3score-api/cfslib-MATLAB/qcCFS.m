function [status, quality, message] = qcCFS(filename, threshold)

data = readCFS(filename);
quality = sum(squeeze(mean(mean(squeeze(data))))>800,2)*100/size(data,4);
status = 0;
electrodes = {'C3/C4', 'EoG-L', 'EoG-R', 'EMG'};
message = 'All channels passed quality checks.';
failed_channels = '';
qc = quality > threshold;
if(any(qc))
    status = 1;
    for i = 1: numel(qc),
        if qc(i) 
            failed_channels = strcat(failed_channels, electrodes{i}, ', ');
        end
    end
    message = strcat('The following channel(s) failed quality checks: ', failed_channels);
end

end