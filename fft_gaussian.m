close all;
clear;

%generate the time index
steps = 1000;
dt = 1;%tmax / steps;
tmax = dt*steps;
t = dt:dt:tmax;

%determine the frequency of the input signal
F = 1000;

Tm = dt * 5;
T0 = tmax / 2;%Tm * 10;
x = exp(-(t - T0).^2/Tm/dt*1E-2); %.*sin(2*F*pi*t);

figure;
plot(t, x);


ff = 1:steps;
%apply the FFT transform on the input signals
y = fftshift(fft(x(ff).*hamming(length(ff))'));
%  y = fft(x);

%generate the frequency index
f = (ff - (max(ff) + min(ff))/2 ) / (length(ff) *dt);

%plot the frequency components of the input signal
figure;
plot(f, abs(y));
xlabel('Frequency (Hz)'); 
ylabel('Abs. Magnitude'); grid on;

figure;
plot(f, (abs(real(y))));

return;

% % % % % % % % % % % % % % % % % % % % % % % % 
%Continuous Wavelet Transform
figure;
c = cwt(x, 0.01:1:100, 'db4', 'plot');
colormap(jet);

% % % % % % % % % % % % % % % % % % % % % % % % % % % % %%%%
% sum c along time should be the Fourier frequency spectrum 
figure;
plot(abs(sum(c')));
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% calculate the average value of c along t
avg_steps = 10;
c_avg = zeros(size(c,1), round(steps/avg_steps));
for i = 1:avg_steps:steps
    
end;


% % % % % % % % % % % % % % % % % % % % % % % % % % % 
% Short-time Fourier Transform of x
figure;
Fs=1E3;                  % sampled every 0.001 sec so rate is 1 kHz
step=ceil(10*Fs/1000);    % one spectral slice every 20 ms
window=ceil(100*Fs/1000); % 100 ms data window
%specgram(x, 2^nextpow2(window), Fs, window, window-step);
specgram(x, 2^nextpow2(window), Fs, window, window-step);

% % % % % % % % % % % % % % % % % % % % % % % % % % % 
% Short-time Fourier Transform 
figure;
xx = chirp([0:0.001:2],0,2,500);  % freq. sweep from 0-500 over 2 sec.
Fs=1000;                  % sampled every 0.001 sec so rate is 1 kHz
step=ceil(10*Fs/1000);    % one spectral slice every 20 ms
window=ceil(100*Fs/1000); % 100 ms data window
specgram(xx, 2^nextpow2(window), Fs, window, window-step);

