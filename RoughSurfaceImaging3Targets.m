%% RoughSurfaceImaging3Targets.m
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
% generate Figs. 10 and 11 in that manuscript. 
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

% target 1

x1   = -9.0;
z1   = -10.1;
rho1 = 3.6 * 1j;

% target 2

x2   =  1.0;
z2   = -9.4;
rho2 = 3.4 * 1j;

% target 3

x3   = 11.0;
z3   = -9.8;
rho3 = 3.6 * 1j;

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

        Rsca1 = sqrt( ( x - x1 ).^2 + ( h - z1 ).^2 );
        Usca1 = 1j / 4 * besselh( 0, 1, k1 * Rsca1 );

        Rsca2 = sqrt( ( x - x2 ).^2 + ( h - z2 ).^2 );
        Usca2 = 1j / 4 * besselh( 0, 1, k1 * Rsca2 );

        Rsca3 = sqrt( ( x - x3 ).^2 + ( h - z3 ).^2 );
        Usca3 = 1j / 4 * besselh( 0, 1, k1 * Rsca3 );

        % solve the system of BIEs

        c = [ A11 A21; A21 A22 ] \ ...
            [ Usrc.' zeros(N,1) zeros(N,1) zeros(N,1); zeros(N,1) Usca1.' Usca2.' Usca3.' ];
        
        % parse the solutions
        
        U0   = c(1:N,1);
        DnU0 = c(N+1:2*N,1);

        U1   = c(1:N,2);
        DnU1 = c(N+1:2*N,2);

        U2   = c(1:N,3);
        DnU2 = c(N+1:2*N,3);

        U3   = c(1:N,4);
        DnU3 = c(N+1:2*N,4);

        % compute the reflected field on z = za

        Ra   = sqrt( ( xa(j) - x ).^2 + ( za - h ).^2 );           
        SLPa = 1j / 4 * besselh( 0, k0 * Ra );
        DLPa = 1j * k0 / 4 * ( za - h ) .* besselh( 1, k0 * Ra ) ./ Ra;

        Dref(i,j) = dx * ( DLPa * U0 - SLPa * DnU0 );
        
        % compute the transmitted field on z = z1 and z2

        Rtarget1   = sqrt( ( x1 - x ).^2 + ( z1 - h ).^2 );
        SLPtarget1 = 1j / 4 * besselh( 0, k1 * Rtarget1 );
        DLPtarget1 = 1j * k1 / 4 * ( z1 - h ) ...
            .* besselh( 1, k1 * Rtarget1 ) ./ Rtarget1;
        
        Rtarget2   = sqrt( ( x2 - x ).^2 + ( z2 - h ).^2 );
        SLPtarget2 = 1j / 4 * besselh( 0, k1 * Rtarget2 );
        DLPtarget2 = 1j * k1 / 4 * ( z2 - h ) ...
            .* besselh( 1, k1 * Rtarget2 ) ./ Rtarget2;

        Rtarget3   = sqrt( ( x3 - x ).^2 + ( z3 - h ).^2 );
        SLPtarget3 = 1j / 4 * besselh( 0, k1 * Rtarget3 );
        DLPtarget3 = 1j * k1 / 4 * ( z3 - h ) ...
            .* besselh( 1, k1 * Rtarget3 ) ./ Rtarget3;
        
        Uexcite1 = - dx * ( DLPtarget1 * U0 - SLPtarget1 * DnU0 );
        Uexcite2 = - dx * ( DLPtarget2 * U0 - SLPtarget2 * DnU0 );
        Uexcite3 = - dx * ( DLPtarget3 * U0 - SLPtarget3 * DnU0 );
        
        % compute the scattered field on z = za

        Dsca(i,j) = rho1 * Uexcite1 * dx * ( DLPa * U1 - SLPa * DnU1 ) ...
            + rho2 * Uexcite2 * dx * ( DLPa * U2 - SLPa * DnU2 ) ...
            + rho3 * Uexcite3 * dx * ( DLPa * U3 - SLPa * DnU3 );

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

% remove the first principal component from the data

Dtilde = DD - sigma(1) * U(:,1) * V(:,1)' ...
    - sigma(2) * U(:,2) * V(:,2)' ...
    - sigma(3) * U(:,3) * V(:,3)' ...
    - sigma(4) * U(:,4) * V(:,4)' ...
    - sigma(5) * U(:,5) * V(:,5)';

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
plot( x1, z1, 'r.', 'MarkerSize', 6 );
plot( x1, z1, 'ro', 'MarkerSize', 14, 'LineWidth', 1 );
plot( x2, z2, 'r.', 'MarkerSize', 6 );
plot( x2, z2, 'ro', 'MarkerSize', 14, 'LineWidth', 1 );
plot( x3, z3, 'r.', 'MarkerSize', 6 );
plot( x3, z3, 'ro', 'MarkerSize', 14, 'LineWidth', 1 );
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
plot( x1, z1, 'r.', 'MarkerSize', 6 );
plot( x1, z1, 'ro', 'MarkerSize', 14, 'LineWidth', 1 );
plot( x2, z2, 'r.', 'MarkerSize', 6 );
plot( x2, z2, 'ro', 'MarkerSize', 14, 'LineWidth', 1 );
plot( x3, z3, 'r.', 'MarkerSize', 6 );
plot( x3, z3, 'ro', 'MarkerSize', 14, 'LineWidth', 1 );
hold off;
xlabel( '$x$ [cm]', 'Interpreter', 'LaTeX' );
ylabel( '$z$ [cm]', 'Interpreter', 'LaTeX' );

%% LOCAL IMAGE ABOUT TARGET 1

KM1  = KM(abs(zIW-z1)<=5,abs(xIW-x1)<=5);
rKM1 = epsilon ./ ( 1 - (1 - epsilon) * KM1 / max(KM1(:)) );
xIW1 = xIW(abs(xIW-x1)<=5);
zIW1 = zIW(abs(zIW-z1)<=5);

figure(9)
pcolor( kcentral * ( xIW1 - x1 ), kcentral * ( zIW1 - z1 ), rKM1 );
shading interp;
colorbar;
clim( [ 0 1 ] );
axis equal tight;
hold on;
plot( 0, 0, 'r.', 'MarkerSize', 6 );
plot( 0, 0, 'ro', 'MarkerSize', 14, 'LineWidth', 1 );
hold off;
xlabel( '$k_{0}( x - x_{1} )$', 'Interpreter', 'LaTeX' );
ylabel( '$k_{0}( z - z_{1} )$', 'Interpreter', 'LaTeX' );

%% LOCAL IMAGE ABOUT TARGET 2

KM2  = KM(abs(zIW-z2)<=5,abs(xIW-x2)<=5);
rKM2 = epsilon ./ ( 1 - (1 - epsilon) * KM2 / max(KM2(:)) );
xIW2 = xIW(abs(xIW-x2)<=5);
zIW2 = zIW(abs(zIW-z2)<=5);

figure(10)
pcolor( kcentral * ( xIW2 - x2 ), kcentral * ( zIW2 - z2 ), rKM2 );
shading interp;
colorbar;
clim( [ 0 1 ] );
axis equal tight;
hold on;
plot( 0, 0, 'r.', 'MarkerSize', 6 );
plot( 0, 0, 'ro', 'MarkerSize', 14, 'LineWidth', 1 );
hold off;
xlabel( '$k_{0}( x - x_{2} )$', 'Interpreter', 'LaTeX' );
ylabel( '$k_{0}( z - z_{2} )$', 'Interpreter', 'LaTeX' );

%% LOCAL IMAGE ABOUT TARGET 3

KM3  = KM(abs(zIW-z3)<=5,abs(xIW-x3)<=5);
rKM3 = epsilon ./ ( 1 - (1 - epsilon) * KM3 / max(KM3(:)) );
xIW3 = xIW(abs(xIW-x3)<=5);
zIW3 = zIW(abs(zIW-z3)<=5);

figure(11)
pcolor( kcentral * ( xIW3 - x3 ), kcentral * ( zIW3 - z3 ), rKM3 );
shading interp;
colorbar;
clim( [ 0 1 ] );
axis equal tight;
hold on;
plot( 0, 0, 'r.', 'MarkerSize', 6 );
plot( 0, 0, 'ro', 'MarkerSize', 14, 'LineWidth', 1 );
hold off;
xlabel( '$k_{0}( x - x_{3} )$', 'Interpreter', 'LaTeX' );
ylabel( '$k_{0}( z - z_{3} )$', 'Interpreter', 'LaTeX' );