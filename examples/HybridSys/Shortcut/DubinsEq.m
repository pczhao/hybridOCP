function xdot = DubinsEq( t, x, controller, J )

polysin = @(x) x;
polycos = @(x) 1 - x^2/2;

xval = x(1:3);
uval = controller( t, xval );
uval(uval>1) = 1;
uval(uval<0) = 0;

if nargin < 4
    xdot = [ uval(1) * polycos(xval(3));
             uval(1) * polysin(xval(3));
             3 * uval(2) ];
else
    xdot = [ uval(1) * polycos(xval(3));
             uval(1) * polysin(xval(3));
             3 * uval(2);
             J(xval) ];
end