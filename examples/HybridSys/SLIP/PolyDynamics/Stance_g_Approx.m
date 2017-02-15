function [out] = Stance_g_Approx( ~, params )
k = params.k;
m = params.m;
umax = params.umax;

out = [ 0;
        k/m * umax;
        0;
        0;
        0 ];