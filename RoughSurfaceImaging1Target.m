%% RoughSurfaceImaging1Target.m
%
% Solve the system of boundary integral equations for adjacent half spaces
% separated by a rough interface to compute measurements for synthetic 
% aperture subsurface imaging. Then apply PCA and Kirchhoff migration to 
% solve the inverse scattering problem. This code uses the Fourier 
% quadrature method to compute illuminations for the KM illuminations.
%
% This is the code used to generate the images in the manuscript,
% "Synthetic aperture radar imaging below a random rough surface," by A. D.
% Kim and C. Tsogka. The parameters in the code below are those used to 
% generate Figs. 4, 5, and 6 in that manuscript. The left plot in Fig. 6 is
% obtained by replacing "Dtilde" with "Dsca" in line 371 below.
%
% Written by A. D. Kim on 3/20/2023

clear;

%% FIGURE PARAMETERS

set(0,'defaultaxesfontsize',20,'defaultaxeslinewidth',1.0,...
      'defaultlinelinewidth',2.0,'defaultpatchlinewidth',1.0); 

%% PHYSICAL PARAMETERS

% wavenumbers

Nf   = 25;
f    = linspace( 3.1, 5.1, Nf );  % 1 GHz - 4 GHz
c0   = 29.9792458;                % speed of light in centimeters GHz
k    = 2.0 * pi * f / c0;         % wavenumbers
nrel = 3.00;                      % relative refractive index
beta = 0.10;                      % loss tangent

kcentral = mean( k );

% mean interface location

zbdy = 0.0;

% source locations

Na = 21;
a  = 100;
xa = linspace( - a / 2, a / 2, Na );
za = 100;

% target location and reflectivity

xtarget =  2.0;
ztarget = -8.0;
rho     = 3.4 * 1j;

%% X GRID

N  = 512;
L  = 400;
dx = L / N;
x  = - L / 2 : dx : L / 2 - dx;
xk = 2.0 * pi / L * fftshift( - N / 2 : N / 2 - 1 );

%% ROUGH SURFACE

rng( 'default' );

hRMS = 0.20;
lcor = 8.00;

h   = zbdy + GenerateRandomFunction( N, L, hRMS, lcor );
dh  = ifft( -1j * xk .* fft( h ) );
d2h = ifft( -xk.^2 .* fft( h ) );

%% PLOT THE REALIZATION OF THE ROUGH SURFACE

figure(1)
plot( kcentral * x, kcentral * h );
axis( [ -100 100 -2 2 ] );
grid on;
xlabel( '$k_{0} x$', 'Interpreter', 'LaTeX' );
ylabel( '$k_{0} h(x)$', 'Interpreter', 'LaTeX' );

%% COMPUTE MEASUREMENTS

% allocate memory for the measurements

Dref = zeros( Nf, Na );
Dsca = zeros( Nf, Na );

% define useful index arrays

[ indx, jndx ] = ndgrid( (1:N) );

% compute useful matrices for SLP and DLP operators

R = sqrt( ( x(indx) - x(jndx) ).^2 + ( h(indx) - h(jndx) ).^2 );
V = dh(jndx) .* ( x(indx) - x(jndx) ) - ( h(indx) - h(jndx) );

% approximate Euler-Mascheroni Constant

gamma = 1.78107;

% loop over frequencies

disp( ' ' );
disp( 'Computing measurements' );
disp( ' ' );

tic

for i = 1 : Nf

    % set wavenumbers

    k0 = k(i);
    k1 = k(i) * nrel * sqrt( 1 + 1j * beta );
    
    % compute BIE operators
    
    S0 = 0 * R;
    S1 = 0 * R;
    D0 = 0 * R;
    D1 = 0 * R;
    
    % non-diagonal entries

    S0(R~=0) = dx * 1j / 4 * besselh( 0, k0 * R(R~=0) );
    S1(R~=0) = dx * 1j / 4 * besselh( 0, k1 * R(R~=0) );
    D0(R~=0) = dx * 1j * k0 / 4 * besselh( 1, k0 * R(R~=0) ) ./ R(R~=0) ...
        .* V(R~=0);
    D1(R~=0) = dx * 1j * k1 / 4 * besselh( 1, k1 * R(R~=0) ) ./ R(R~=0) ...
        .* V(R~=0);

    % diagonal entries

    S0 = S0 + diag( dx * 1j / 4 * ( 1 + 1j * 2 / pi ...
            * log( gamma * k0 / 4 / exp(1) * dx * sqrt( 1 + dh.^2 ) ) ) );

    S1 = S1 + diag( dx * 1j / 4 * ( 1 + 1j * 2 / pi ...
            * log( gamma * k1 / 4 / exp(1) * dx * sqrt( 1 + dh.^2 ) ) ) );

    D0 = D0 + diag( -d2h / ( 4 * pi ) * dx ./ ( 1 + dh.^2 ) );

    D1 = D1 + diag( -d2h / ( 4 * pi ) * dx ./ ( 1 + dh.^2 ) );

    % construct system of BIEs

    A11 = 0.5 * eye(N) - D0;
    A12 = S0;
    A21 = 0.5 * eye(N) + D1;
    A22 = -S1;

    % loop over spatial locations

    for j = 1 : Na

        % compute sources

        Rsrc = sqrt( ( x - xa(j) ).^2 + ( h - za ).^2 );
        Usrc = 1j / 4 * besselh( 0, 1, k0 * Rsrc );

        Rsca = sqrt( ( x - xtarget ).^2 + ( h - ztarget ).^2 );
        Usca = 1j / 4 * besselh( 0, 1, k1 * Rsca );

        % solve the system of BIEs

        c = [ A11 A21; A21 A22 ] \ [ Usrc.' zeros(N,1); zeros(N,1) Usca.' ];
        
        % parse the solutions
        
        U0   = c(1:N,1);
        DnU0 = c(N+1:2*N,1);

        U1   = c(1:N,2);
        DnU1 = c(N+1:2*N,2);

        % compute the reflected field on z = za

        Ra   = sqrt( ( xa(j) - x ).^2 + ( za - h ).^2 );           
        SLPa = 1j / 4 * besselh( 0, k0 * Ra );
        DLPa = 1j * k0 / 4 * ( za - h ) .* besselh( 1, k0 * Ra ) ./ Ra;

        Dref(i,j) = dx * ( DLPa * U0 - SLPa * DnU0 );

        % compute the transmitted field on z = ztarget

        Rtarget   = sqrt( ( xtarget - x ).^2 + ( ztarget - h ).^2 );
        SLPtarget = 1j / 4 * besselh( 0, k1 * Rtarget );
        DLPtarget = 1j * k1 / 4 * ( ztarget - h ) ...
            .* besselh( 1, k1 * Rtarget ) ./ Rtarget;
        
        Uexcite = - dx * ( DLPtarget * U0 - SLPtarget * DnU0 );
        
        % compute the scattered field on z = za

        Dsca(i,j) = rho * Uexcite * dx * ( DLPa * U1 - SLPa * DnU1 );

    end

end

toc

%% MEASUREMENT NOISE

sn     = 1e-2;
D      = Dref + Dsca;
D0     = D;
enr    = sqrt(sum(abs(D(:)).^2)/length(D(:)));
Dnoise = enr * sn * (randn(size(D)) + 1j * randn(size(D))) / sqrt(2);
DD     = D0 + Dnoise;
 
SNR_est  = -10 * log10( sn );
SNR_real =  10 * log10( norm(D0) / norm(Dnoise) );
SNR_sca  =  10 * log10( norm(Dsca) / norm(Dnoise) );

disp( ' ' );
disp( '---' );
disp( 'SNR' );
disp( '---' );
disp( ' ' );
disp( ['  Realization SNR = ', num2str(SNR_real), ' dB' ] );
disp( ['  Average SNR     = ', num2str(SNR_est), ' dB'  ] );
disp( ['  Effective SNR   = ', num2str(SNR_sca), ' dB'  ] );
disp( ' ' );

%% GROUND BOUNCE REMOVAL USING PCA

% compute the SVD

[ U, Sigma, V ] = svd( DD );

% plot the singular values

sigma = diag( Sigma );

figure(2)

semilogy( (1:length(sigma)), sigma / sigma(1), '.', 'MarkerSize', 12 );
grid on;
ylim( [ 1e-4 10 ] );
xlabel( '$j$', 'Interpreter', 'LaTeX' );
ylabel( '$\sigma_{j}/\sigma_{1}$', 'Interpreter', 'LaTeX' );

% remove the first rtrunc principal components from the data

Dtilde = DD;

rtrunc = 5;

for j = 1 : rtrunc

    Dtilde = Dtilde - sigma(j) * U(:,j) * V(:,j)';

end

%% plot the processed measurements

% measurements

figure(3)

pcolor( xa, f, real(DD) );
shading flat;
colorbar;
xlabel( '$x_{n}$ [cm]', 'Interpreter', 'LaTeX' );
ylabel( '$f$ [GHz]', 'Interpreter', 'LaTeX' );
% title( 'raw measurements', 'Interpreter', 'LaTeX' );

% ground bounce

figure(4)

pcolor( xa, f, real(Dref) );
shading flat;
colorbar;
xlabel( '$x_{n}$ [cm]', 'Interpreter', 'LaTeX' );
ylabel( '$f$ [GHz]', 'Interpreter', 'LaTeX' );
% title( 'ground bounce', 'Interpreter', 'LaTeX' );

% scattered

figure(5)

pcolor( xa, f, real(Dsca) );
shading flat;
colorbar;
xlabel( '$x_{n}$ [cm]', 'Interpreter', 'LaTeX' );
ylabel( '$f$ [GHz]', 'Interpreter', 'LaTeX' );
% title( 'scattered', 'Interpreter', 'LaTeX' );

% PCA processed measurements

figure(6)

pcolor( xa, f, real(Dtilde) );
shading flat;
colorbar;
xlabel( '$x_{n}$ [cm]', 'Interpreter', 'LaTeX' );
ylabel( '$f$ [GHz]', 'Interpreter', 'LaTeX' );
% title( 'processed', 'Interpreter', 'LaTeX' );

%% COMPUTE ILLUMINATIONS

% xz grid for the imaging window

NN     = 101;
xIW    = linspace( -15, 15, NN );
zIW    = linspace( -20, -1, NN );

[ X, Z ] = meshgrid( xIW, zIW );

X = X(:);
Z = Z(:);

% set the order of the Fourier quadrature rule

M = 500;

% useful index arrays

[ indx, jndx ] = ndgrid( 1 : NN * NN, 1 : M );

% allocate memory for KM imaging function

KM = 0 * X;

% message to the user

disp( ' ' );
disp( 'Computing KM imaging function' );
disp( ' ' );

% loop over frequencies and spatial locations

tic

for i = 1 : Nf

    % wavenumbers

    k0 = k(i);
    k1 = k(i) * nrel;
    
    % Fourier quadrature rule

    [ xi, wt ] = FourierQuadrature( k0, k1, za, M );

    % propagation constants

    q0 = sqrt( k0^2 - xi(jndx).^2 );
    q1 = sqrt( k1^2 - xi(jndx).^2 );

    for j = 1 : Na

        % field from source to target

        U1hat = 0.5 * 1j / pi ./ ( q0 + q1 ) ...
            .* exp( 1j * q0 * za - 1j * q1 .* Z(indx) ) ...
            .* exp( 1j * xi(jndx) .* ( X(indx) - xa(j) ) );

        u1 = U1hat * wt.';

        % field from target to source

        U0hat = 0.5 * 1j / pi ./ ( q0 + q1 ) ...
            .* exp( 1j * q0 * za - 1j * q1 .* Z(indx) ) ...
            .* exp( -1j * xi(jndx) .* ( X(indx) - xa(j) ) );

        u0 = U0hat * wt.';

        % update KM imaging function

        KM = KM + Dtilde(i,j) * conj( u0 ./ abs( u0 ) ) ...
            .* conj( u1 ./ abs( u1 ) ); 

    end

end

toc

%% KM

KM = reshape( abs( KM ), NN, NN );
KM = KM / max( KM(:) );

% rKM

epsilon = 1e-2;
rKM     = epsilon ./ ( 1 - ( 1 - epsilon ) * KM );

%% PLOT THE IMAGES

figure(7)
pcolor( xIW, zIW, KM );
shading flat;
colorbar;
clim( [ 0 1 ] );
axis equal tight;
hold on;
plot( xtarget, ztarget, 'r.', 'MarkerSize', 6 );
plot( xtarget, ztarget, 'ro', 'MarkerSize', 14, 'LineWidth', 1 );
hold off;
xlabel( '$x$ [cm]', 'Interpreter', 'LaTeX' );
ylabel( '$z$ [cm]', 'Interpreter', 'LaTeX' );

figure(8)
pcolor( xIW, zIW, rKM );
shading flat;
colorbar;
clim( [ 0 1 ] );
axis equal tight;
hold on;
plot( xtarget, ztarget, 'r.', 'MarkerSize', 6 );
plot( xtarget, ztarget, 'ro', 'MarkerSize', 14, 'LineWidth', 1 );
hold off;
xlabel( '$x$ [cm]', 'Interpreter', 'LaTeX' );
ylabel( '$z$ [cm]', 'Interpreter', 'LaTeX' );
