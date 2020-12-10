% function to create CFS stream (byte array) from raw EEG data
% raw data is 4 X Sample matrix, the 4 channels in order are:
% C3, C4, EL and ER respectively
%
% (c)-2018 Neurobit Technologies - Amiya Patanaik 
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
% THE SOFTWARE.
function EDF2CFS(filepath)
    if(nargin ~= 1),
        [filename,path,~] = uigetfile('*.edf');
        filepath = fullfile(path,filename);
    end
    [path,filename,~] = fileparts(filepath);
    [header, signalHeader, signalCell] = blockEdfLoad([path '/' filename '.edf']);
    N = numel(signalHeader);
    disp('Here are all the channels:');
    for k=1:N,
        fprintf('%d: %s\n',k,signalHeader(k).signal_labels);
    end
    C3N = input('Please select the C3:A2 channel number: ');
    C4N = input('Please select the C4:A1 channel number: ');
    ELN = input('Please select the EOGl:A2 channel number: ');
    ERN = input('Please select the EOGr:A1 channel number: ');
    EMGN = input('Please select the chin EMG channel number: ');
    
    f_EEG = signalHeader(C3N).samples_in_record/header.data_record_duration;
    f_EOG = signalHeader(ELN).samples_in_record/header.data_record_duration;
    f_EMG = signalHeader(EMGN).samples_in_record/header.data_record_duration;
    
    samplingRate = signalHeader(C3N).samples_in_record/header.data_record_duration;
    EEGData = [signalCell{C3N}'; signalCell{C4N}'; signalCell{ELN}'; signalCell{ERN}'];
    
    writeCFS([path '/' filename '_v1.cfs'], EEGData, samplingRate);
    writeCFS_V2([path '/' filename '_v2.cfs'], signalCell{C3N}, signalCell{C4N}, signalCell{ELN}, signalCell{ERN}, signalCell{EMGN}, f_EEG, f_EOG, f_EMG);
    
end
    
