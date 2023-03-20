function [ f ] = GenerateRandomFunction( N, rL, h, lc )
%
% ------------------------------
% MATLAB: GenerateRoughSurface.m
% ------------------------------
%
% Generate a realization of a 1D Gaussian random rough surfaces 
% with Gaussian Spectrum.
%
% INPUT:
%
% N     = total number of sample points
% rL    = rough surface length
% h     = rms height
% lc    = correlation length
% iseed = seed of random number generator
%
% OUTPUT:
%
% f  = rough surface profile
%

y = randn(N,1);

for n = 1:(N/2-1);
    
  bh(n) = ( y(2*n-1) + i*y(2*n) ) / sqrt(2.0);
  
end;

bhc = conj(bh);
bhf = fliplr(bhc);
bi  = [bh y(N-1) bhf y(N)];
kx  = 2*pi*[-N/2+1:1:N/2]/rL;
y1  = sqrt( wk(kx,h,lc) );
y   = y1*sqrt(2*pi*rL);
b   = y .* bi;
xs  = [ b(N/2+1:1:N) b(1:1:N/2) ];
xt  = [xs(N),xs(1:1:N-1)];
ft  = ifft(xt,N);
ft  = ft * N / rL;
fs  = [ft(2:1:N),ft(1)];
f   = [fs(N/2+1:1:N) fs(1:1:N/2)];
f   = real(f);

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gaussian spectral density %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = wk(kx,h,lc)

y=h^2*lc*exp(-(kx*lc*0.5).^2)/(2*sqrt(pi));
