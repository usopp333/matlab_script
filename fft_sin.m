close all;
clear;

%generate the time index
tmax = 2;
dt = 0.0001;
steps = tmax / dt;
t = dt:dt:tmax;

fprintf(1, 'steps:%d\n', steps);
%determine the frequency of the input signal
F = 200;

x = sin(2*F*pi*t);

alpha = 1E1;
x2 = sin(2*F*pi*t).*exp(-alpha*t);
x2_padding = x2;
x2_padding(length(x2)/10:length(x2)) = 0;
x2_trunc = x2(1:length(x2)/10);

figure;
%plot(t, x, t, x2, t, x2_padding);
plot(t, x2, t, x2_padding);
legend('No Padding', 'With Padding');

figure;
%plot(t, x, t, x2, t, x2_padding);
plot(t, x2);


ff = 1:steps;
%apply the FFT transform on the input signals
y = fftshift(fft(x(ff)));
y2 = fftshift(fft(x2(ff)));
y2_padding = fftshift(fft(x2_padding));
y2_trunc = fftshift(fft(x2_trunc));

y3 = fftshift(fft(x(ff).*hamming(length(ff))'));
y4 = fftshift(fft(x2(ff).*hamming(length(ff))'));
%  y = fft(x);

%generate the frequency index
f = (ff - (max(ff) + min(ff))/2 ) / (length(ff) *dt);

ff2 = 1:length(x2)/10;
f_trunc = (ff2 - (max(ff2) + min(ff2))/2 ) / (length(ff2) *dt);


y_analytical = alpha/(2*pi) ./ ( (f - F).^2  + (alpha/(2*pi))^2 );



%plot the frequency components of the input signal
figure;
%plot(f, abs(y), f, abs(y2), f, abs(y3), f, abs(y4));
%legend('sin', 'lossy sin', 'sin with hamming window',  'lossy sin with hamming window');
plot(f, abs(y2), f, abs(y4));
%legend('lossy sin',  'lossy sin with hamming window');
legend('Without hamming window',  'With hamming window');
xlabel('Frequency (Hz)'); 
ylabel('Abs. Magnitude'); grid on;




ratio = max(abs(y2)) / max(y_analytical);

figure;
plot(f, abs(y2), f, abs(y2_padding), f_trunc, abs(y2_trunc));
legend('No Padding', 'With Padding', 'Truncated');
xlabel('Frequency (Hz)'); 
ylabel('Abs. Magnitude'); grid on;




% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
figure_fontsize = 12;
set(get(gca,'Xlabel'), 'FontSize', figure_fontsize, 'Vertical', 'top');
set(get(gca,'Ylabel'), 'FontSize', figure_fontsize, 'Vertical', 'bottom');
%set(findobj('FontSize', 10), 'FontSize', figure_fontsize);
set(findobj(get(gca, 'Children'), 'Linewidth', 0.5), 'Linewidth', 2);




