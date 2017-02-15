function [ R ] = Reset_F2S_Approx( x, params )
% we are not approximating cos(alpha) or sin(alpha).
alpha = params.alpha;
l0 = params.l0;

polysin = @(xx) sin(xx);
polycos = @(xx) cos(xx);

R = [ l0;
      -x(2) * polysin(alpha) + x(4) * polycos(alpha);
      alpha;
      -x(2) * polycos(alpha) / l0 - x(4) * polysin(alpha) / l0;
      x(1) ];
