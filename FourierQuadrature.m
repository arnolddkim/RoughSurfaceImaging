%% FourierQuadrature.m
%
% This function computes the quadrature rule on a deformed contour for
% computing Fourier integrals.
%
% Written by A. D. Kim on 12/3/2021

function [ xi, wt ] = FourierQuadrature( k0, k1, zmin, M )

% evaluate a conformal map introduced by Barnett and Greengard to cluster 
% quadrature points in the interval [-k0,k0]

lambda = 1.0;
delta  = 1.0;
K      = sqrt( ( log( eps ) / zmin )^2 + k0^2 );
L      = lambda * asinh( K / ( delta * lambda ) );

% compute the beta grid

beta = -L + 2 * L * ( 1 : M ) / M;

% conformal map parameter

s = lambda * delta * sinh( beta / lambda );

% deformed contour

A = 0.40;
w = 6.00;

if k0 == k1
    
    A = A * 0.5;
    
end

% quadrature points

xi = s + 1j * A * ( exp( -w * ( s + k0 ).^2 ) ...
    + exp( -w * ( s + k1 ).^2 ) ...
    - exp( -w * ( s - k0 ).^2 ) ...
    - exp( -w * ( s - k1 ).^2 ) );

% quadrature weights

wt = 1.0 - 1j * 2 * A * w * ( ( s + k0 ) .* exp( -w * ( s + k0 ).^2 ) ...
    + ( s + k1 ) .* exp( -w * ( s + k1 ).^2 ) ...
    - ( s - k0 ) .* exp( -w * ( s - k0 ).^2 ) ...
    - ( s - k1 ) .* exp( -w * ( s - k1 ).^2 ) );

wt = wt * 2 * L / M * delta .* cosh( beta / lambda );

return;