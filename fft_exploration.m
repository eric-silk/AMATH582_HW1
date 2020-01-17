close all; clear all; clc
n = 256;
L = 1;

t = linspace(-L, L, n);
freq = 8;
y = sin(2*pi*freq*t);

figure();
subplot(311);
plot(t,y);

k = (1/(2*L))*[0:(n/2-1) -n/2:-1];
ks = fftshift(k);
Y = fftshift(fft(y));
subplot(312);
stem(ks, abs(Y))

yt = ifft(fftshift(Y));
subplot(313)
plot(t, yt)

figure();
hist(abs(Y))
