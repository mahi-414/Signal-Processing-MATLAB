% Lab 7: Applications of Laplace Filtering, Transfer Functions, FFTs

%{
In this lab, we applied what we learned about the above topics to convert
an MRI brain image from the frequency domain to the spatial domain, as well
as clean up a noisy sound signal.
%}

%% Filtering MRI Brain Image
clear all;

% loading data
load("EBME358_Lab7_MRI_data.mat")

% processing data
imageData = abs(ifft2(fftshift(MRIdata)));

% plotting original data
figure
imagesc(imageData);
colormap gray;
axis off;
axis equal;
title("MRI image")

% plotting original data on subplot
figure
subplot(2,2,1)
imagesc(imageData);
colormap gray;
axis off;
axis equal;
title("No cols/rows zero")

% plotting original data with every other column zeroed
MRIdata2 = MRIdata;
MRIdata2(:,1:2:end) = 0;
imageData2 = abs(ifft2(fftshift(MRIdata2)));
subplot(2,2,2)
imagesc(imageData2);
colormap gray;
axis off;
axis equal;
title("Every other col=0")

% plotting original data with every other row zeroed
MRIdata3 = MRIdata;
MRIdata3(1:2:end,:) = 0;
imageData3 = abs(ifft2(fftshift(MRIdata3)));
subplot(2,2,3)
imagesc(imageData3);
colormap gray;
axis off;
axis equal;
title("Every other row=0")

% plotting original data with every other row and column zeroed
MRIdata4 = MRIdata2;
MRIdata4(1:2:end,:) = 0;
imageData4 = abs(ifft2(fftshift(MRIdata4)));
subplot(2,2,4)
imagesc(imageData4);
colormap gray;
axis off;
axis equal;
title("Every other row and col = 0")

% low pass filter radius 50
low_pass50 = zeros(128);
low_pass50(15:114, 15:114) = 1;
MRIdata_lp50 = MRIdata .* low_pass50;
imageData_lp50 = abs(ifft2(fftshift(MRIdata_lp50)));

figure
imagesc(imageData_lp50);
colormap gray;
axis off;
axis equal;
title("low pass filter = 50")

% low pass filter radius 25
low_pass25 = zeros(128);
low_pass25(39:88, 39:88) = 1;
MRIdata_lp25 = MRIdata .* low_pass25;
imageData_lp25 = abs(ifft2(fftshift(MRIdata_lp25)));

figure
imagesc(imageData_lp25);
colormap gray;
axis off;
axis equal;
title("low pass filter = 25")

% high pass filter radius 50
high_pass50 = ones(128);
high_pass50(15:114, 15:114) = 0;
MRIdata_hp50 = MRIdata .* high_pass50;
imageData_hp50 = abs(ifft2(fftshift(MRIdata_hp50)));

figure
subplot(2,1,1)
imagesc(imageData_hp50);
colormap gray;
axis off;
axis equal;
title("high pass filter = 50")

% high pass filter radius 25
high_pass25 = ones(128);
high_pass25(39:88, 39:88) = 0;
MRIdata_hp25 = MRIdata .* high_pass25;
imageData_hp25 = abs(ifft2(fftshift(MRIdata_hp25)));

subplot(2,1,2)
imagesc(imageData_hp25);
colormap gray;
axis off;
axis equal;
title("high pass filter = 25")


%% Filtering a real-world signal
% Fs = 44100 Hz
p = load('EBME358_Lab7_P4.mat');
signal = p.Q2;

figure

% signal in time domain
Fs = 44100;
deltaT = 1/Fs;
T0 = (length(signal)-1)*deltaT;
t = 0:deltaT:T0;

subplot(3,1,1)
plot(t, signal)
title('Signal in Time Domain')
xlabel('time (s)')
ylabel('signal amplitude')

% signal in freq domain
freq = -Fs/2:1/T0:Fs/2;
DFT = abs(fftshift(fft(signal*deltaT)));

subplot(3,1,2)
plot(freq, DFT)
title('Signal in Freq Deomain')
xlabel('frequency (Hz)')
ylabel('signal amplitude')

% spectrogram
subplot(3,1,3)
window = 1000;
overlap = 500;
f = 0:Fs/2; 
spectrogram(signal, window, overlap, f, Fs, 'yaxis')
title("Spectrogram")

colorbar off
colormap bone

%% LPF and HPF of order 1, Wc = 5000 and 500

% LOW PASS FILTERS
figure
Wc = 5000; % 5000 Hz
[num, denom] = butter(1, Wc, 'low', 's');
LPF = tf(num, denom);
subplot(3,2,1)
bp = bodeplot(LPF);
opts = getoptions(bp);
opts.FreqUnits='Hz';
opts.MagUnits='abs';
opts.Grid='on';
setoptions(bp,opts);
legend("1st Order")
title('Bode plot of Butterworth LPF Wc = 5000Hz')

[n, d] = bilinear(num, denom, Fs);
LPF_signal = filter(n,d, signal);
subplot(3,2,3)
plot(t, LPF_signal)
title('Signal after Butterworth LPF Wc = 5000Hz')

freq = -Fs/2:1/T0:Fs/2;
DFT = abs(fftshift(fft(LPF_signal*deltaT)));

subplot(3,2,5)
plot(freq, DFT)
title('Signal after Butterworth LPF in Freq Deomain')
xlabel('frequency (Hz)')
ylabel('signal amplitude')


% HIGH PASS FILTERS
Wc = 500; % 500 Hz
[num, denom] = butter(1, Wc, 'high', 's');
HPF = tf(num, denom);
subplot(3,2,2)
bp = bodeplot(HPF);
opts = getoptions(bp);
opts.FreqUnits='Hz';
opts.MagUnits='abs';
opts.Grid='on';
setoptions(bp,opts);
legend("1st Order")
title('Bode plot of Butterworth HPF Wc = 500Hz')

[n, d] = bilinear(num, denom, Fs);
HPF_signal = filter(n,d, LPF_signal);
subplot(3,2,4)
plot(t, HPF_signal)
title('Signal after Butterworth HPF Wc = 500Hz')

freq = -Fs/2:1/T0:Fs/2;
DFT = abs(fftshift(fft(HPF_signal*deltaT)));

subplot(3,2,6)
plot(freq, DFT)
title('Signal after Butterworth HPF in Freq Deomain')
xlabel('frequency (Hz)')
ylabel('signal amplitude')

%% LPF and HPF of order 1, Wc = 3000 and 1000

% LOW PASS FILTERS
figure
Wc = 3000; % 3000 Hz
[num, denom] = butter(1, Wc, 'low', 's');
LPF = tf(num, denom);
subplot(3,2,1)
bp = bodeplot(LPF);
opts = getoptions(bp);
opts.FreqUnits='Hz';
opts.MagUnits='abs';
opts.Grid='on';
setoptions(bp,opts);
legend("1st Order")
title('Bode plot of Butterworth LPF Wc = 3000Hz')

[n, d] = bilinear(num, denom, Fs);
LPF_signal = filter(n,d, signal);
subplot(3,2,3)
plot(t, LPF_signal)
title('Signal after Butterworth LPF Wc = 3000Hz')

freq = -Fs/2:1/T0:Fs/2;
DFT = abs(fftshift(fft(LPF_signal*deltaT)));

subplot(3,2,5)
plot(freq, DFT)
title('Signal after Butterworth LPF in Freq Deomain')
xlabel('frequency (Hz)')
ylabel('signal amplitude')


% HIGH PASS FILTERS
Wc = 1000; % 1000 Hz
[num, denom] = butter(1, Wc, 'high', 's');
HPF = tf(num, denom);
subplot(3,2,2)
bp = bodeplot(HPF);
opts = getoptions(bp);
opts.FreqUnits='Hz';
opts.MagUnits='abs';
opts.Grid='on';
setoptions(bp,opts);
legend("1st Order")
title('Bode plot of Butterworth HPF Wc = 1000Hz')

[n, d] = bilinear(num, denom, Fs);
HPF_signal = filter(n,d, LPF_signal);
subplot(3,2,4)
plot(t, HPF_signal)
title('Signal after Butterworth HPF Wc = 1000Hz')

freq = -Fs/2:1/T0:Fs/2;
DFT = abs(fftshift(fft(HPF_signal*deltaT)));

subplot(3,2,6)
plot(freq, DFT)
title('Signal after Butterworth HPF in Freq Deomain')
xlabel('frequency (Hz)')
ylabel('signal amplitude')

%% LPF and HPF of order 5, Wc = 3000 and 1000

% LOW PASS FILTERS
figure
Wc = 3000; % 3000 Hz
[num, denom] = butter(5, Wc, 'low', 's');
LPF = tf(num, denom);
subplot(3,2,1)
bp = bodeplot(LPF);
opts = getoptions(bp);
opts.FreqUnits='Hz';
opts.MagUnits='abs';
opts.Grid='on';
setoptions(bp,opts);
legend("5th Order")
title('Bode plot of Butterworth LPF Wc = 3000Hz')

[n, d] = bilinear(num, denom, Fs);
LPF_signal = filter(n,d, signal);
subplot(3,2,3)
plot(t, LPF_signal)
title('Signal after Butterworth LPF Wc = 3000Hz')

freq = -Fs/2:1/T0:Fs/2;
DFT = abs(fftshift(fft(LPF_signal*deltaT)));

subplot(3,2,5)
plot(freq, DFT)
title('Signal after Butterworth LPF in Freq Deomain')
xlabel('frequency (Hz)')
ylabel('signal amplitude')


% HIGH PASS FILTERS
Wc = 1000; % 1000 Hz
[num, denom] = butter(5, Wc, 'high', 's');
HPF = tf(num, denom);
subplot(3,2,2)
bp = bodeplot(HPF);
opts = getoptions(bp);
opts.FreqUnits='Hz';
opts.MagUnits='abs';
opts.Grid='on';
setoptions(bp,opts);
legend("5th Order")
title('Bode plot of Butterworth HPF Wc = 1000Hz')

[n, d] = bilinear(num, denom, Fs);
HPF_signal = filter(n,d, LPF_signal);
subplot(3,2,4)
plot(t, HPF_signal)
title('Signal after Butterworth HPF Wc = 1000Hz')

freq = -Fs/2:1/T0:Fs/2;
DFT = abs(fftshift(fft(HPF_signal*deltaT)));

subplot(3,2,6)
plot(freq, DFT)
title('Signal after Butterworth HPF in Freq Deomain')
xlabel('frequency (Hz)')
ylabel('signal amplitude')

% FINAL SIGNAL
final_signal = HPF_signal;

%%
sound(final_signal, Fs)