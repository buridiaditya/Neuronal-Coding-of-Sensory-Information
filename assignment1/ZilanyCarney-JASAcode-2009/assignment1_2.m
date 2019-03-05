
%
CF = 500; % CF in Hz;
cohc  = 1.0;   % normal ohc function
cihc  = 1.0;   % normal ihc function
fiberType = 1; % spontaneous rate (in spikes/s) of the fiber BEFORE refractory effects; "1" = Low; "2" = Medium; "3" = High
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
intensityRange = -20:10:200;

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
experimentData = ones(1,length(intensityRange));

figure;
hold on;
% ANF model 
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Fs any arbitrary value does not work
Fs = 100e3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;
for j=1:length(intensityRange)
    intensity = intensityRange(j);
    pin = clippedAudio*sqrt(2)*20e-6*10^(intensity/20.0)/Rmsvalue;
    [synout, psth] = ANModel(nrep, pin, CF, Fs, length(pin)/Fs, cohc, cihc, fiberType,implnt); 
    experimentData(1,j) = sum(psth);
end


plot(intensityRange, experimentData,'DisplayName','ah vowel tone')

for j=1:length(intensityRange)
    intensity = intensityRange(j);
    pin = generateStimulus(CF, Fs, length(clippedAudio)/Fs, rt, intensity);
    [synout, psth] = ANModel(nrep, pin, CF, Fs, length(pin)/Fs, cohc, cihc, fiberType,implnt); 
    experimentData(1,j) = sum(psth);
end
toc;
plot(intensityRange, experimentData,'DisplayName','BF tone')

legend()
%%
% 40, 100, 190

sound1 = sqrt(2)*20e-6*10^(40/20.0)* audio /Rmsvalue;
sound2 = sqrt(2)*20e-6*10^(100/20.0)* audio /Rmsvalue;
sound3 = sqrt(2)*20e-6*10^(190/20.0)* audio /Rmsvalue;

%%
window = 25.6e-3*Fs;
noverlap = window/2;
figure;

spectrogram(audio, window, noverlap, window, Fs, 'yaxis')
set(gca,'Yscale','log')

%%
figure;

windows = [4,8,16,32,64,128];
tic;
for index=1:length(windows)
    window = windows(index)*1e-3*Fs;
    overlap = window/2;

    [synout, psth] = ANModel(nrep, sound1', frequencyBag(1), Fs, length(sound1)/Fs, cohc, cihc, fiberType,implnt);     
    psth = processPSTH(psth, window, overlap);
    experimentData = zeros(length(frequencyBag), length(psth));
    experimentData(1,:) = psth;

    
    
    % ANF model 1
    parfor i=2:length(frequencyBag)
        CF = frequencyBag(i);
        [synout, psth] = ANModel(nrep, sound1', CF, Fs, length(sound1)/Fs, cohc, cihc, fiberType,implnt);     
        psth = processPSTH(psth, window, overlap);
        experimentData(i,:) = psth;
    end
    

    subplot(3,2,index) 
    
    imagesc(gca, 'XData', [1,length(psth)], 'YData', [frequencyBag(1), frequencyBag(length(frequencyBag))], 'CData', flipud(experimentData))
    title( strcat(int2str(window),'ms windowsize'))
    ytick = 80*2.^(0:1/8:7);
    set(gca,'YTick',ytick);

end
toc;
%%
figure;

windows = [4,8,16,32,64,128];
tic;
for index=1:length(windows)
    window = windows(index)*1e-3*Fs;
    overlap = window/2;

    [synout, psth] = ANModel(nrep, sound2', frequencyBag(1), Fs, length(sound2)/Fs, cohc, cihc, fiberType,implnt);     
    psth = processPSTH(psth, window, overlap);
    experimentData = zeros(length(frequencyBag), length(psth));
    experimentData(1,:) = psth;

    
    
    % ANF model 1
    parfor i=2:length(frequencyBag)
        CF = frequencyBag(i);
        [synout, psth] = ANModel(nrep, sound2', CF, Fs, length(sound2)/Fs, cohc, cihc, fiberType,implnt);     
        psth = processPSTH(psth, window, overlap);
        experimentData(i,:) = psth;
    end
    

    subplot(3,2,index) 
    
    imagesc(gca, 'XData', [1,length(psth)], 'YData', [frequencyBag(1), frequencyBag(length(frequencyBag))], 'CData', flipud(experimentData))
    title( strcat(int2str(window),'ms windowsize'))
    ytick = 80*2.^(0:1/8:7);
    set(gca,'YTick',ytick);

end
toc;

%%
figure;

windows = [4,8,16,32,64,128];
tic;
for index=1:length(windows)
    window = windows(index)*1e-3*Fs;
    overlap = window/2;

    [synout, psth] = ANModel(nrep, sound3', frequencyBag(1), Fs, length(sound3)/Fs, cohc, cihc, fiberType,implnt);     
    psth = processPSTH(psth, window, overlap);
    experimentData = zeros(length(frequencyBag), length(psth));
    experimentData(1,:) = psth;

    
    
    % ANF model 1
    parfor i=2:length(frequencyBag)
        CF = frequencyBag(i);
        [synout, psth] = ANModel(nrep, sound3', CF, Fs, length(sound3)/Fs, cohc, cihc, fiberType,implnt);     
        psth = processPSTH(psth, window, overlap);
        experimentData(i,:) = psth;
    end
        
    subplot(3,2,index) 
    
    imagesc(gca, 'XData', [1,length(psth)], 'YData', [frequencyBag(1), frequencyBag(length(frequencyBag))], 'CData', flipud(experimentData))
    title( strcat(int2str(window),'ms windowsize'))
    ytick = 80*2.^(0:1/8:7);
    set(gca,'YTick',ytick);
end
toc;
