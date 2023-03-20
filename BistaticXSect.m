%% BistaticXsect.m
%
% Compute the backscattering cross-section for the rough surfaces used in 
% the subsurface imaging code for the manuscript entitled, "Synthetic
% aperture radar imaging below a random rough surface" by A. D. Kim and C.
% Tsogka. The code is set to produce the image used in the left plot of
% Fig. 2 of that manuscript.
%
% This code implements the method described in Chapter 4 of Tsang, Kong, 
% Ding, and Ao, "Scattering of Electromagnetic Waves: Numerical 
% Simulations" (Wiley, 2001).
%
% Written by A. D. Kim on 3/20/2023

clear;

%% FIGURE PARAMETERS

set(0,'defaultaxesfontsize',20,'defaultaxeslinewidth',1.0,...
      'defaultlinelinewidth',2.0,'defaultpatchlinewidth',1.0); 

%% PHYSICAL PARAMETERS

% wavenumbers

Nf   = 25;
f    = linspace( 3.1, 5.1, Nf );  % 3.1 GHz - 5.1 GHz
c0   = 29.9792458;                % speed of light in centimeters GHz
k    = 2.0 * pi * f / c0;         % wavenumbers
nrel = 3.00;                      % relative refractive index
beta = 0.10;                      % loss tangent

% mean interface location

zbdy = 0.0;

%% X GRID

N  = 512;
L  = 400;
dx = L / N;
x  = - L / 2 : dx : L / 2 - dx;
xk = 2.0 * pi / L * fftshift( - N / 2 : N / 2 - 1 );

%% ROUGH SURFACE PARAMETERS

hRMS = 0.20;
lcor = 8.00;

%% INCIDENT FIELD -- TAPERED PLANE WAVE PARAMETERS

g       = L / 4;
Nang    = 101;
theta_i = 30 * pi / 180;
theta_s = theta_i + linspace( -25, 25, Nang ) * pi / 180;

%% COMPUTE BISTATIC CROSS-SECTION

% set gamma constant

gamma = 1.78107;

% allocate memory for bistatic cross-section

sigma = zeros( Nf, Nang );

% define useful index arrays

[ indx, jndx ] = ndgrid( (1:N) );
[ iang, jx   ] = ndgrid( (1:Nang), (1:N) );

% seed random number generator

rng( 'default' );

% number of realizations

Nsample = 100;

% loop over realizations

for isamp = 1 : Nsample

    % generate rough surface realization
    
    h   = zbdy + GenerateRandomFunction( N, L, hRMS, lcor );
    dh  = ifft( -1j * xk .* fft( h ) );
    d2h = ifft( -xk.^2 .* fft( h ) );
    
    % compute useful matrices for SLP and DLP operators
    
    R  = sqrt( ( x(indx) - x(jndx) ).^2 + ( h(indx) - h(jndx) ).^2 );
    V = dh(jndx) .* ( x(indx) - x(jndx) ) - ( h(indx) - h(jndx) );
    
    % loop over frequencies
    
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
    
        % incident field -- tapered plane wave
    
        w = ( 2 * ( x + h * tan( theta_i ) ).^2 / g^2  - 1 ) ...
            ./ ( k0 * g * cos( theta_i ) ).^2;
        
        a = exp( - ( x + h * tan( theta_i ) ).^2 / g^2 );
    
        b1 = exp( 1j * k0 * ( x * sin( theta_i ) - h * cos( theta_i ) ) ...
            .* ( 1 + w ) ) .* a;
    
        b2 = 0 * b1;
    
        % solve the system of BIEs
    
        c = [ A11 A21; A21 A22 ] \ [ b1.'; b2.' ];
            
        % parse the solutions
        
        U   = c(1:N,1);
        DnU = c(N+1:2*N,1);
    
        % scattering amplitude
    
        psiN = dx * sum( -DnU(jx) + U(jx) * 1j * k0 ...
            .* ( dh(jx) .* sin( theta_s(iang) ) - cos( theta_s(iang) ) ) ...
            .* exp( -1j * k0 * ( x(jx) .* sin( theta_s(iang) ) ...
            + h(jx) .* cos( theta_s(iang) ) ) ), 2 );
    
        % bistatic cross-section
    
        sigma(i,:) = sigma(i,:) + abs( psiN' ).^2 ...
            ./ ( 8 * pi * k0 * g * sqrt( pi / 2 ) ...
            * cos( theta_i ) * ( 1 - ( 1 + 2 * tan( theta_i )^2 ) ...
            / ( 2 * k0^2 * g^2 * cos( theta_i )^2 ) ) ) / Nsample;
    
    end

    disp( [ '   Computed realization ', num2str(isamp), ' out of ', num2str(Nsample) ] );

end

%% PLOT THE RESULTS

mesh( theta_s - theta_i, f, 10 * log10( sigma ) )
xlabel( '$\theta_{s}$ [rad]', 'Interpreter', 'LaTeX' );
ylabel( '$f$ [GHz]', 'Interpreter', 'LaTeX' );
zlabel( '$\langle \sigma(\theta_{s}) \rangle$ [dB]', 'Interpreter', 'LaTeX' );
colormap( [ 0 0 0 ] );
view( -18, 41 );