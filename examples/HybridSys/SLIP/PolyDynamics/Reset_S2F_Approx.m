function [ out ] = Reset_S2F_Approx( x, params )
l0 = params.l0;

polysin = @(x) x;
polycos = @(x) 1 - x^2/2;

out = [ x(5);
        -x(2) * polysin(x(3)) - l0 * x(4) * polycos(x(3));
        l0 * polycos(x(3));
        x(2) * polycos(x(3)) - l0 * x(4) * polysin(x(3)) ];
