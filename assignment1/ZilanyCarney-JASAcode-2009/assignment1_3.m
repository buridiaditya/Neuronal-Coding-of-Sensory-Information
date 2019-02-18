
%
CF = 500; % CF in Hz;
cohc  = 1.0;   % normal ohc function
cihc  = 1.0;   % normal ihc function
fiberType = 2; % spontaneous rate (in spikes/s) of the fiber BEFORE refractory effects; "1" = Low; "2" = Medium; "3" = High
implnt = 0;    % "0" for approximate or "1" for actual implementation of the power-law functions in the Synapse
%

%
% stimulus parameters
rt = 10e-3;   % rise/fall time in seconds

% PSTH parameters
nrep = 50;               % number of stimulus repetitions (e.g., 50);
psthbinwidth = 0.5e-3; % binwidth in seconds;

%
frequencyBag = 80*2.^(0:1/8:7);
intensityRange = -50:10:200;
%%

[audio, Fs] = audioread('fivewo.wav');
audioinfo('fivewo.wav')

T = length(audio)*(1/Fs);

clippedAudio = audio(110000:1:120000)';

%for i=1:100
%    sound(clippedAudio, Fs);
%end
 
Rmsvalue = rms(clippedAudio);
dBspl = 20 * log10(Rmsvalue/20*10^(-6));

    
mxpts = length(clippedAudio);
irpts = double(int64(rt*Fs));

clippedAudio(1:irpts) = clippedAudio(1:irpts).*(0:(irpts-1))/irpts; 
clippedAudio((mxpts-irpts):mxpts) = clippedAudio((mxpts-irpts):mxpts).*(irpts:-1:0)/irpts;


%%
Fs = 100e3;
window = 12.8e-3*Fs;
overlap = window/2;

[synout, psth1] = ANModel(nrep, audio', CF, Fs, length(audio)/Fs, cohc, cihc, fiberType,implnt);     
psth1 = processPSTH(psth1, window, overlap);
experimentData = zeros(length(frequencyBag), length(psth1));

%%
tic;
% ANF model 1
lenF = length(frequencyBag);
parfor i=1:length(frequencyBag)
    CF = frequencyBag(i);
    [synout, psth] = ANModel(nrep, audio', CF, Fs, length(audio1)/Fs, cohc, cihc, fiberType,implnt);     
    psth = processPSTH(psth, window, overlap);
    experimentData(i,:) = psth;
end
toc;
figure;
imagesc('XData', [1,length(psth1)], 'YData', [frequencyBag(1), frequencyBag(length(frequencyBag))], 'CData', flipud(experimentData))

%% 
figure;
hold on;
plot(audio);
for i=1:length(frequencyBag)
    transform = fft(audio);    
end
    