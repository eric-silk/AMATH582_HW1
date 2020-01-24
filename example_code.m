clear all; close all; clc;

% Configuration Constants
% Note that plotting all isosurfaces will take considerable
% time -- 5 minutes is not unrealistic.
% If PLOT_ALL_*_SURFACES is set to false, the user must
% press enter to continue after each plot
PLOT_TIME_ISOSURFACE = false;
PLOT_ALL_TIME_SURFACES = false;
PLOT_FREQ_ISOSURFACE = false;
PLOT_ALL_FREQ_SURFACES = false;

load Testdata
SAMPLES = 20;
L=15; % spatial domain
n=64; % Fourier modes
x2=linspace(-L,L,n+1);
x=x2(1:n); y=x; z=x;
k=(2*pi/(2*L))*[0:(n/2-1) -n/2:-1];
ks=fftshift(k); % Seems like an odd way but ok
[X,Y,Z]=meshgrid(x,y,z);
[Kx,Ky,Kz]=meshgrid(ks,ks,ks);

% Time Domain
if PLOT_TIME_ISOSURFACE
    for j=1:SAMPLES
        figure();
        Un(:,:,:)=reshape(Undata(j,:),n,n,n);
        isosurface(X,Y,Z,abs(Un),1.25)
        axis([-20 20 -20 20 -20 20]), grid on
        if not PLOT_ALL_TIME_SURFACES
            drawnow;
            input("Press ENTER to continue...");
            close all;
        end
    end
    drawnow;
end

% Frequency Domain
if PLOT_FREQ_ISOSURFACE
    for j=1:SAMPLES
        figure();
        Un(:,:,:)=reshape(Undata(j,:),n,n,n);
        Un_fft=fftn(Un);
        isosurface(Kx,Ky,Kz,abs(fftshift(Un_fft)))
        axis([-2*pi 2*pi -2*pi 2*pi -2*pi 2*pi])
        grid on
        if not PLOT_ALL_FREQ_SURFACES
            drawnow;
            input("Press ENTER to contiue...");
            close all;
        end
    end
    drawnow;
end

averaged_spectrum = zeros(n, n, n);
figure
% Plot histograms of the frequency magnitudes while summing
% all the spectrums
% This will inform our cutoff for the isosurface
for j=1:SAMPLES
    tmp(:,:,:) = reshape(Undata(j, :), n, n, n);
    fft_tmp = fftn(tmp);
    subplot(4,5,j)
    % Note: not norming yet!
    histogram(abs(fft_tmp))
    title(sprintf("Sample %d", j))
    averaged_spectrum = averaged_spectrum + fft_tmp;
end
sgtitle({"Spectrum Histograms", "X axis: Absolute Magnitude, Y axis: Count"})

% Now get the average instead of just sum
averaged_spectrum = abs(averaged_spectrum) / SAMPLES;
figure
histogram(averaged_spectrum)
title("Averaged Spectrum")
xlabel("Spectrum Magnitude")
ylabel("Count")

plot_spectrum_isosurface(Kx, Ky, Kz, averaged_spectrum, 100);
title("Averaged Spectrum, Surface=100")
plot_spectrum_isosurface(Kx, Ky, Kz, averaged_spectrum, 125);
title("Averaged Spectrum, Surface=125")
plot_spectrum_isosurface(Kx, Ky, Kz, averaged_spectrum, 150);
title("Averaged Spectrum, Surface=150")
plot_spectrum_isosurface(Kx, Ky, Kz, averaged_spectrum, 175);
title("Averaged Spectrum, Surface=175")
plot_spectrum_isosurface(Kx, Ky, Kz, averaged_spectrum, 200);
title("Averaged Spectrum, Surface=200")
plot_spectrum_isosurface(Kx, Ky, Kz, averaged_spectrum, 250);
title("Averaged Spectrum, Surface=250")

% Generate a 3D gaussian centered about x=1.885, y=-1.08, z=0
sigma_x = 1;
sigma_y = 1;
sigma_z = 1;
mu_x = 1.885;
mu_y = -1.08;
mu_z = 0;

coeff = 1/(2*pi*sigma_x*sigma_y*sigma_z);
exponent = (Kx-mu_x).^2/(2*sigma_x^2) + (Ky-mu_y).^2/(2*sigma_y^2) + (Kz-mu_z).^2/(2*sigma_z^2);
gauss_filter = coeff*exp(-exponent);
% normalize it
gauss_filter = gauss_filter/(max(gauss_filter(:)));
figure();
title(sprintf("Filter kernel\nsigma=%.1f, mu_x=%.3f, mu_y=%.3f, mu_z=%.3f", sigma_x, mu_x, mu_y, mu_z));
isosurface(Kx, Ky, Kz, gauss_filter, 0.9);
axis([-2*pi 2*pi -2*pi 2*pi -2*pi 2*pi]);
xlabel("Frequency, X")
ylabel("Frequency, Y")
zlabel("Frequency, Z")
grid on, drawnow

figure();
for j=1:20
    tmp(:,:,:) = reshape(Undata(j, :), n, n, n);
    fft_tmp = fftshift(fftn(tmp));
    filtered = fft_tmp.*gauss_filter;
    filtered_time = ifftn(fftshift(filtered));
    subplot(4,5,j);
    histogram(abs(filtered_time));
    title(sprintf("Sample %d", j))
end
sgtitle({"Histograms of Filtered Spectrums",...
         "X axis: Absolute Magnitude, Y axis: Count"})

% Colors, dark to light for time
figure();
title("Progress of Marble in Fluffy's innards")
xlabel("X")
ylabel("Y")
zlabel("Z")
hold on;

num = SAMPLES;
den = num + 1;
init_color = num / den;
color_iter = -[1/den 1/den 1/den];
facecolor = [init_color init_color init_color];
for j=1:SAMPLES
    tmp(:,:,:) = reshape(Undata(j, :), n, n, n);
    fft_tmp = fftshift(fftn(tmp));
    filtered = fft_tmp .* gauss_filter;
    filtered_time = ifftn(fftshift(filtered));
    hiso = patch(isosurface(X, Y, Z, abs(filtered_time), 0.3),...
                 'FaceColor', facecolor,...
                 'EdgeColor', 'none');
    isonormals(abs(filtered_time), hiso);
    axis([-20 20 -20 20 -20 20]);
    grid on;
    facecolor = facecolor + color_iter;
end

% Now lets do the plot3 trajectory plot
x3d = zeros(SAMPLES);
y3d = zeros(SAMPLES);
z3d = zeros(SAMPLES);
for j=1:SAMPLES
    tmp(:,:,:) = reshape(Undata(j, :), n, n, n);
    fft_tmp = fftshift(fftn(tmp));
    filtered = fft_tmp .* gauss_filter;
    filtered_time = ifftn(fftshift(filtered));
    [max_val, max_index] = max(filtered_time(:));
    [indx, indy, indz] = ind2sub(size(filtered_time), max_index);
    if not (max(tmp(:)) == tmp(indx, indy, indz))
        disp("Ya done goofed, A-A-ron!")
    end
    x3d(j) = indx;
    y3d(j) = indy;
    z3d(j) = indz;
end

figure
plot3(x3d, y3d, z3d)

disp("It's over...it's done.")