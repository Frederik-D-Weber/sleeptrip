% function to create CFS version 2 stream (byte array) from raw PSG data
% CFS V2 supports the NEO sleep classification system, which has shown 35
% to 40% lower error in many datasets as compared to Z3Score
% the 5 channels in order are: C3-A2, C4-A1, EL-A2, ER-A1 and EMG respectively
% 
% streamCFS_V2(EEGData, EOGData, EMG, samplingRate, varargin);
% EEGData is raw 2 channel raw EEG data comprising of C3-A2, C4-A1. 
% fsamp_EEG is the sampling rate, must be sampled at 100Hz or more
%
% EOGData is raw 2 channel raw EOG data comprising of EL-A2, ER-A1. 
% fsamp_EOG is the sampling rate, must be sampled at 100Hz or more
%
% EMGData is raw 1 channel raw EMG data comprising of chin EMG. 
% fsamp_EMG is the sampling rate, must be sampled at 200Hz or more
%
% Compress is optional argument (default = 1), to control compression
% Hash is optional argument (default = 1), to control transport security
%
% (c)-2018 Neurobit Technologies Amiya Patanaik amiya@neurobit.io
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
% THE SOFTWARE.

function [stream, status, message, quality] = streamCFS_V2(C3, C4, EOGL, EOGR, EMG, fsamp_EEG, fsamp_EOG, fsamp_EMG, varargin)

% set defaults for optional inputs
optargs = {1 1};
numvarargs = length(varargin);
if numvarargs > 8,
    error('streamCFS accepts 10 arguments only.')
end
optargs(1:numvarargs) = varargin;
[compress, hash] = optargs{:};

% Check sampling rates
if fsamp_EEG < 100
    disp('Sampling rate must be > 100Hz for EEG channels');
    return
end

if fsamp_EOG < 100
    disp('Sampling rate must be > 100Hz for EOG channels');
    return
end

if fsamp_EMG < 200
    disp('Sampling rate must be > 200Hz for EMG channels');
    return
end

%Order 50 FIR filter
%Basic Settings
LOWPASS = 35.0;  % Hz
HIGHPASS = 0.3;  % Hz
LOWPASSEOG = 35.0;  % Hz
LOWPASSEMG = 80.0; % Hz
threshold = 10;

%DC blocking filter EEG
p = dcblock(0.1,fsamp_EEG);             
b = [1 -1];                         % set up differentiator
a = [1 -p];                         % set up integrator
C3 = filter(b,a,C3); 
C4 = filter(b,a,C4);

%DC blocking filter EOG
p = dcblock(0.1,fsamp_EOG);             
b = [1 -1];                         % set up differentiator
a = [1 -p];                         % set up integrator
EOGL = filter(b,a,EOGL); 
EOGR = filter(b,a,EOGR);

%DC blocking filter EMG
p = dcblock(0.1,fsamp_EMG);             
b = [1 -1];                         % set up differentiator
a = [1 -p];                         % set up integrator
EMG = filter(b,a,EMG); 



Fs_EEG = fsamp_EEG/2;
Fs_EOG = fsamp_EOG/2;
Fs_EMG = fsamp_EMG/2;


bEEG = fir1(50,[HIGHPASS/Fs_EEG LOWPASS/Fs_EEG]);
bEOG = fir1(50,[HIGHPASS/Fs_EOG LOWPASSEOG/Fs_EOG]);
bEMG = fir1(50,[HIGHPASS/Fs_EMG LOWPASSEMG/Fs_EMG]);

eogL = filter(bEOG,1,EOGL);
eogR = filter(bEOG,1,EOGR);
eeg = (filter(bEEG,1,C3) + filter(bEEG,1,C4))./2;
emg = filter(bEMG,1,EMG);

if fsamp_EEG ~= 100,
    [p,q] = rat(100/fsamp_EEG);
    eeg =  resample(eeg,p,q);
end

if fsamp_EOG ~= 100,
    [p,q] = rat(100/fsamp_EOG);
    eogL = resample(eogL,p,q);
    eogR = resample(eogR,p,q);
end

if fsamp_EMG ~= 200,
    [p,q] = rat(200/fsamp_EMG);
    emg =  resample(emg,p,q);
end

%Extract features from data in each epoch

totalEpochs = floor(size(eeg,1)/30/100);
data = zeros(32,32,4,totalEpochs);

EPOCH = 30*100; %Samples
EPOCH_EMG = 30*200; %Samples

fun = @(block_struct) mean(block_struct.data(:));

for j=1:totalEpochs,
    
    % FOR EEG-------------------------------------------
    s = spectrogram(eeg((j-1)*EPOCH+1:j*EPOCH),128,36,128,100,'yaxis');
    data(:,:,1,j) = abs(s(2:33,:));
       
    % FOR EOGL-------------------------------------------
    s = spectrogram(eogL((j-1)*EPOCH+1:j*EPOCH),128,36,128,100,'yaxis');
    data(:,:,2,j) = abs(s(2:33,:));
    
    % FOR EOGR-------------------------------------------
    s = spectrogram(eogR((j-1)*EPOCH+1:j*EPOCH),128,36,128,100,'yaxis');
    data(:,:,3,j) = abs(s(2:33,:));
    
    % FOR EMG-------------------------------------------
    s = spectrogram(emg((j-1)*EPOCH_EMG+1:j*EPOCH_EMG),256,71,256,200,'yaxis');
    s = abs(s(2:end,:));
    data(:,:,4,j) = blockproc(s, [4 1], fun);
    
end

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

%Signature first 3 bytes
signature = uint8('CFS');
%file version next 1 byte
version = uint8(2);
%frequency 1 byte - time 1byte - channel 1byte epochs 2bytes
frequency = uint8(32);
time = uint8(32);
channel = uint8(4);
epochs = typecast(uint16(totalEpochs),'uint8');
%boolean bits
compressionbit = uint8(compress);
hashsetbit = uint8(hash);

%Convert data to single as per specifications
data = cast(data(:),'single')';
actualDigest = [];

if(hashsetbit == 1)
    Opt.Method = 'SHA-1';
    Opt.Format = 'uint8';
    Opt.Input = 'bin';
    actualDigest = DataHash(data, Opt);
end

if(compressionbit == 1)
    data = zlibencode(typecast(data,'uint8'));   
else
    data = typecast(data,'uint8'); 
end

stream = [signature version frequency time channel epochs compressionbit hashsetbit actualDigest data];

end