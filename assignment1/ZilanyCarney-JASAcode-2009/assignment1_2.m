
%
CF = 500; % CF in Hz;
CF1 = 4000;
cohc  = 1.0;   % normal ohc function
cihc  = 1.0;   % normal ihc function
fiberType = 3; % spontaneous rate (in spikes/s) of the fiber BEFORE refractory effects; "1" = Low; "2" = Medium; "3" = High
implnt = 0;    % "0" for approximate or "1" for actual implementation of the power-law functions in the Synapse
%

%
% stimulus parameters
F0 = CF;     % stimulus frequency in Hz
Fs = 10e3;  % sampling rate in Hz (must be 100, 200 or 500 kHz)
T  = 2;  % stimulus duration in seconds
rt = 10e-3;   % rise/fall time in seconds

% PSTH parameters
nrep = 50;               % number of stimulus repetitions (e.g., 50);
psthbinwidth = 0.5e-3; % binwidth in seconds;

%

frequencyBag = 62.5*2.^[0:1/8:7];

[audio, Fs] = audioread('fivewo.wav');

% Find the segment "ah"
% How to get "ah" ?
clippedAudio ;

Rmsvalue = rms(clippedAudio);

dBspl = 20 * log10(Rmsvalue/20*10^(-6));

intensityRange = [-150:10:150];

experimentData = ones(1,length(intensityRange))

figure
% ANF model 1
for j=1:length(intensityRange)
    intensity = intensityRange(j);
    pin = clippedAudio*sqrt(2)*20e-6*10^(intensity/20.0)/Rmsvalue
    %% How to get T ? 
    [synout, psth] = ANModel(nrep, pin, CF, Fs, T, cohc, cihc, fiberType,implnt); 
    experimentData(1,j) = sum(psth);
end
plot(intensityRange, experimentData)

sound1 = sqrt(2)*20e-6*10^(intensity1/20.0)* audio /Rmsvalue;
sound2 = sqrt(2)*20e-6*10^(intensity2/20.0)* audio /Rmsvalue;
sound3 = sqrt(2)*20e-6*10^(intensity3/20.0)* audio /Rmsvalue;


%%

% ANF model 1
for i=1:length(frequencyBag)
    CF = freqRange(i);
    intensity = intensityRange(j);
    [synout, psth] = ANModel(nrep, sound1, CF, Fs, T, cohc, cihc, fiberType,implnt); 
    plot(psth)    
end

%%%%%%%%

nfft = 4096
win = 25.6*e-3

% Spectrogram
[s,f,t] = stft(audio, win, win/2, nfft, Fs)

plot(t, f, s) %% Figure out how to use color

for i=1:length(frequencyBag)
    CF = freqRange(i);
    intensity = intensityRange(j);
    [synout, psth] = ANModel(nrep, sound1, CF, Fs, T, cohc, cihc, fiberType,implnt); 
    plot(t,frequencyBag, psth) %% Figure out how to use color
        
end

%% 