function PlotFrame( state, params, xv, yv, h_2d )
% state is defined as [ l, ldot, theta, thetadot, x, xdot, y, ydot, mode  ]

l0 = params.l0;
alpha = params.alpha;
Mpos = [ state(5), state(7) ];
if state(end) == 1
    Opos = Mpos + [ state(1)*sin(state(3)), -state(1)*cos(state(3)) ];
else
    Opos = Mpos + [ l0*sin(alpha), -l0*cos(alpha) ];
%     Opos = [nan, nan];
end

% Mass
th = 0:pi/50:2*pi;
xvec = state(5) + 0.02*cos(th);
yvec = state(7) + 0.02*sin(th);
patch(xvec, yvec, 1, 'FaceColor', [0.8,0.8,0.8]);

% Spring
[ xvec, yvec ] = springcoord( Mpos, Opos );
if state(end) == 1
    plot(h_2d, xvec, yvec, 'k');
else
%     plot(h_2d, xvec, yvec, 'color', [0.8,0.8,0.8]);
end

% Velocity
if ~isempty(xv)
%     quiver(h_2d, Mpos(1), Mpos(2), 0.06*xv, 0.06*yv, 'color', 'red', 'MaxHeadSize', 1);
end

