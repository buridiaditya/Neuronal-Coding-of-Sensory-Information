
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
frequencyBag1 = 125*2.^(0:1/2:4);
intensityRange = -50:10:200;
%%

[audio, Fs] = audioread('fivewo.wav');
audioinfo('fivewo.wav')

%%
Fs = 100e3;
window = 0.01e-3*Fs;
overlap = floor(window/2);

[synout, psth] = ANModel(nrep, audio', frequencyBag(1), Fs, length(audio)/Fs, cohc, cihc, fiberType,implnt);     
psth = processPSTH(psth, window, overlap);
experimentData = zeros(length(frequencyBag), length(psth));
experimentData1 = zeros(length(frequencyBag1), length(psth));
experimentData(1,:) = psth;

%%
tic;
% ANF model 1
parfor i=2:length(frequencyBag)
    CF = frequencyBag(i);
    [synout, psth] = ANModel(nrep, audio', CF, Fs, length(audio)/Fs, cohc, cihc, fiberType,implnt);     
    psth = processPSTH(psth, window, overlap);
    experimentData(i,:) = psth;
end
toc;
%%
figure;

imagesc('XData', [1,length(psth)], 'YData', [frequencyBag(1), frequencyBag(length(frequencyBag))], 'CData', flipud(experimentData))
ytick = 80*2.^(0:1/8:7);
set(gca,'YTick',ytick);
hold on;
%%
tic;
% ANF model 1
parfor i=1:length(frequencyBag1)
    CF = frequencyBag1(i);
    [synout, psth] = ANModel(nrep, audio', CF, Fs, length(audio)/Fs, cohc, cihc, fiberType,implnt);     
    psth = processPSTH(psth, window, overlap);
    experimentData1(i,:) = psth;
end
toc;

%% 


fftWindow = 12.8e-3*Fs;
fftoverlap = fftWindow/2;
indices = 1:(fftWindow-fftoverlap):length(experimentData1(1,:));

L = fftWindow;                
T = 1/Fs;             % Sampling period       
t = (0:L-1)*T;

L1 = length(experimentData1(1,:));

for i=1:length(frequencyBag1)
    col = rand(1,3);
    for j=1:length(indices)
        if indices(j)+fftWindow-1 > L1
            Y = fft(experimentData1(i,indices(j):end));   
            L = length(experimentData1(i,indices(j):end));                
            T = 1/Fs;             % Sampling period       
            t = (0:L-1)*T;
        else
            Y = fft(experimentData1(i,indices(j):indices(j)+fftWindow-1));   
        end
        P2 = abs(Y/L);
        P1 = P2(1:L/2+1);
        P1(2:end-1) = 2*P1(2:end-1);
        f = Fs*(0:(L/2))/L;
        [argval, argmax] = max(P1);
        if f(argmax) < 20000 && f(argmax) ~= 0
            if j > 1
                plot((indices(j-1) + indices(j))/2, f(argmax),'*','color',col)
            else
                plot((indices(j))/2, f(argmax),'*','color',col)
            end
        end
        
    end
end

%% 


    