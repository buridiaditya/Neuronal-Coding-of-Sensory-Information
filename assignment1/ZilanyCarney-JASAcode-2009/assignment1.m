% model fiber parameters
clear all;
CF = 500; % CF in Hz;
CF1 = 4000;
cohc  = 1.0;   % normal ohc function
cihc  = 1.0;   % normal ihc function
fiberType = 3; % spontaneous rate (in spikes/s) of the fiber BEFORE refractory effects; "1" = Low; "2" = Medium; "3" = High
implnt = 0;    % "0" for approximate or "1" for actual implementation of the power-law functions in the Synapse

% stimulus parameters
F0 = CF;     % stimulus frequency in Hz
Fs = 10e3;  % sampling rate in Hz (must be 100, 200 or 500 kHz)
T  = 200e-3;  % stimulus duration in seconds
rt = 10e-3;   % rise/fall time in seconds
stimdb = 10; % stimulus intensity in dB SPL

% PSTH parameters
nrep = 50;               % number of stimulus repetitions (e.g., 50);
psthbinwidth = 0.5e-3; % binwidth in seconds;

% Experiments
freqRange = 62.5*2.^[0:1.0/8:9];
intensityRange = [-150:10:150];

%% 
%  Experiment 1

experimentData = zeros(1,length(intensityRange));

figure
% ANF model 1
for i=1:length(freqRange)
    for j=1:length(intensityRange)
        freq = freqRange(i);
        intensity = intensityRange(j);
        pin = generateStimulus(freq, Fs, T, rt, intensity);
        [synout, psth] = ANModel(nrep, pin, CF, Fs, T, cohc, cihc, fiberType,implnt); 
        experimentData(1,j) = sum(psth);
    end
    plot(intensityRange, experimentData)
end


figure

% ANF model 1
for i=1:length(freqRange)
    for j=1:length(intensityRange)
        freq = freqRange(i);
        intensity = intensityRange(j);
        stim = generateStimulus(freq, Fs, T, rt, intensity);
        [synout, psth] = ANModel(nrep, stim, CF1, Fs, T, cohc, cihc, fiberType); 
        experimentData(1,j) = sum(psth);
    end
    plot(intensityRange, experimentData)
endfor


