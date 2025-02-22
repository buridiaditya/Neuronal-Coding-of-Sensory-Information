% model fiber parameters
CF = 500; % CF in Hz;
CF1 = 4000;
cohc  = 1.0;   % normal ohc function
cihc  = 1.0;   % normal ihc function
fiberType = 3; % spontaneous rate (in spikes/s) of the fiber BEFORE refractory effects; "1" = Low; "2" = Medium; "3" = High
implnt = 0;    % "0" for approximate or "1" for actual implementation of the power-law functions in the Synapse

% stimulus parameters
F0 = CF;     % stimulus frequency in Hz
Fs = 50e3;  % sampling rate in Hz (must be 100, 200 or 500 kHz)
T  = 200e-3;  % stimulus duration in seconds
rt = 10e-3;   % rise/fall time in seconds
stimdb = 10; % stimulus intensity in dB SPL

% PSTH parameters
nrep = 50;               % number of stimulus repetitions (e.g., 50);
psthbinwidth = 0.5e-3; % binwidth in seconds;

% Experiments
freqRange = 62.5*2.^(0:(1.0/8):9);
intensityRange = -10:10:80;

%% 
%  Experiment 1
% ANF model 1
tic;
experimentData = zeros(length(freqRange),length(intensityRange));

parfor i=1:length(freqRange)
    experimentDataTemp = zeros(1,length(intensityRange));
    for j=1:length(intensityRange)
        freq = freqRange(i);
        intensity = intensityRange(j);
        pin = generateStimulus(freq, Fs, T, rt, intensity);
        [synout, psth] = ANModel(nrep, pin, CF, Fs, T, cohc, cihc, fiberType,implnt); 
        experimentDataTemp(1,j) = sum(psth);
    end
    experimentData(i,:) = experimentDataTemp;    
end
toc;
%%

figure
hold on;

%for i=1:length(freqRange)
%    plot(intensityRange,experimentData(i,:));
%end
xtick = 62.5*2.^(0:(1.0/8):9);
set(gca,'XTick',xtick)
imagesc(gca, [freqRange(1), freqRange(length(freqRange))], [intensityRange(1) , intensityRange(length(intensityRange))] , experimentData');
%imagesc('XData',freqRange,'YData',intensityRange,'CData',experimentData')


%%
% Model 2
tic;
parfor i=1:length(freqRange)
    experimentDataTemp = zeros(1,length(intensityRange));
    for j=1:length(intensityRange)
        freq = freqRange(i);
        intensity = intensityRange(j);
        pin = generateStimulus(freq, Fs, T, rt, intensity);
        [synout, psth] = ANModel(nrep, pin, CF1, Fs, T, cohc, cihc, fiberType,implnt); 
        experimentDataTemp(1,j) = sum(psth);
    end
    experimentData(i,:) = experimentDataTemp;    
end
toc;
%%

figure
hold on;
%for i=1:length(freqRange)
%    plot(intensityRange,experimentData(i,:));
%end


xtick = 62.5*2.^(0:(1.0/8):9);
set(gca,'XTick',xtick)
imagesc(gca,[freqRange(1), freqRange(length(freqRange))], [intensityRange(1) , intensityRange(length(intensityRange))] , experimentData');
%imagesc('XData',freqRange,'YData',intensityRange,'CData',experimentData')

%%

figure
hold on;

for j=1:length(intensityRange)    
    intensity = intensityRange(j);
    pin = generateStimulus(CF, Fs, T, rt, intensity);
    [synout, psth] = ANModel(nrep, pin, CF, Fs, T, cohc, cihc, fiberType,implnt); 
    experimentData(1,j) = sum(psth);
end

plot(intensityRange,experimentData(1,:),'DisplayName','BF=500Hz');

for j=1:length(intensityRange)    
    intensity = intensityRange(j);
    pin = generateStimulus(CF1, Fs, T, rt, intensity);
    [synout, psth] = ANModel(nrep, pin, CF1, Fs, T, cohc, cihc, fiberType,implnt); 
    experimentData(1,j) = sum(psth);
end

plot(intensityRange,experimentData(1,:),'DisplayName','BF=4kHz');
legend();


