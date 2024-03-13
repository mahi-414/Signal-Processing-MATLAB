% Lab 6: Sampling, FFT, IFFT, Spectral Analysis

%{
In this lab, we used our knowledge of the above topics and applied them to
identify the different frequencies of an EEG, using visualization tools
such as Bode plots and spectrograms
%}

load('ArtificialEEG.mat')

% time plot of EEG
figure
subplot(2,1,1)
fs = 1000;
t0 = ((length(EEG)-1/fs)/fs);
time_vec = 0:1/fs:t0;
plot(time_vec, EEG)
xlabel("time (s)")
ylabel("signal values")
title("EEG signal in time domain")

% frequency plot of EEG
subplot(2,1,2)
fft_eeg = abs(fftshift(fft(EEG)));
f_vec = -fs/2:1/t0:fs/2;
plot(f_vec, fft_eeg);
xlim([0 125]);
xlabel("frequency (Hz)")
ylabel("magnitude")
title("EEG signal in frequency domain")

% frequency plots split up by time according to where the signal is in time
% domain
figure
% 0-100 s
subplot(2,3,1)
eeg_0_100 = EEG(1:find_ind(fs, 100));
eeg_fft_0_100 = abs(fftshift(fft(eeg_0_100)));
f_vec_0_100 = -fs/2:1/100:fs/2;
plot(f_vec_0_100, eeg_fft_0_100);
xlim([0 100]);
xlabel("frequency (Hz)")
ylabel("magnitude")
title("FFT EEG time=[0:100]")

% 100-180 s
subplot(2,3,2)
eeg_100_180 = EEG(find_ind(fs, 100):find_ind(fs, 180));
eeg_fft_100_180 = abs(fftshift(fft(eeg_100_180)));
f_vec_100_180 = -fs/2:1/80:fs/2;
plot(f_vec_100_180, eeg_fft_100_180);
xlim([0 100])
xlabel("frequency (Hz)")
ylabel("magnitude")
title("FFT EEG time=[100:180]")

% 180-300 s
subplot(2,3,3)
eeg_180_300 = EEG(find_ind(fs, 180):find_ind(fs, 300));
eeg_fft_180_300 = abs(fftshift(fft(eeg_180_300)));
f_vec_180_300 = -fs/2:1/120:fs/2;
plot(f_vec_180_300, eeg_fft_180_300);
xlim([0 100])
xlabel("frequency (Hz)")
ylabel("magnitude")
title("FFT EEG time=[180:300]")

% 300-500 s
subplot(2,3,4)
eeg_300_500 = EEG(find_ind(fs, 300):find_ind(fs, 500));
eeg_fft_300_500 = abs(fftshift(fft(eeg_300_500)));
f_vec_300_500 = -fs/2:1/200:fs/2;
plot(f_vec_300_500, eeg_fft_300_500);
xlim([0 100])
xlabel("frequency (Hz)")
ylabel("magnitude")
title("FFT EEG time=[300:500]")

% 500-600 s
subplot(2,3,5)
eeg_500_600 = EEG(find_ind(fs, 500):find_ind(fs, 600));
eeg_fft_500_600 = abs(fftshift(fft(eeg_500_600)));
f_vec_500_600 = -fs/2:1/100:fs/2;
plot(f_vec_500_600, eeg_fft_500_600);
xlim([0 100])
xlabel("frequency (Hz)")
ylabel("magnitude")
title("FFT EEG time=[500:600]")

% 630-860 s
subplot(2,3,6)
eeg_630_860 = EEG(find_ind(fs, 630):find_ind(fs, 860));
eeg_fft_630_860 = abs(fftshift(fft(eeg_630_860)));
f_vec_630_860 = -fs/2:1/230:fs/2;
plot(f_vec_630_860, eeg_fft_630_860);
xlim([0 100])
xlabel("frequency (Hz)")
ylabel("magnitude")
title("FFT EEG time=[630:860]")

% spectrogram
figure
window = 2000;
overlap = 1000;
limits = [0:fs/2];

spectrogram(EEG, window, overlap, limits, fs, 'yaxis')
set(gca,"Yscale", "log")
title('Spectrogram')
h = gca;
h.XTickLabel = string(h.XTick * 60);
xlabel('Time (s)');
colormap bone

%% Functions

% find index given time
function f = find_ind(fs, time)
    f=(time+(1/fs))*fs;
end