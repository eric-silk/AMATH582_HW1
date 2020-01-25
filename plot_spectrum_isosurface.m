function plot_spectrum_isosurface(X, Y, Z, unshifted_fft, mag)
   % User must provide titles!
   figure();
   isosurface(X, Y, Z, fftshift(unshifted_fft), mag);
   axis([min(X(:)) max(X(:)) min(Y(:)) max(Y(:)) min(Z(:)) max(Z(:))]);
   xlabel("X (cycle per distance)")
   ylabel("Y (cycle per distance)")
   zlabel("Z (cycle per distance)")
   grid on
   drawnow
end