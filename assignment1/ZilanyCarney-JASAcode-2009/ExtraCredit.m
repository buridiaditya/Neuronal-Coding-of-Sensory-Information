[audio, Fs] = audioread('fivewo.wav');
audioinfo('fivewo.wav')

T = length(audio)*(1/Fs);
%%