function xdot = DIEq( t, x, controller )

uval = controller( t, x );
xdot = [ x(2);
         uval ];
end
