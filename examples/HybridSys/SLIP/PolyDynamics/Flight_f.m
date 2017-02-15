function [out] = Flight_f( x, params )
g = params.g;

% Vectorize
if size(x,1) == 1
    x = x';
end
out = [ x(2,:);
        0 * x(2,:);
        x(4,:);
        -g + 0 * x(2,:) ];