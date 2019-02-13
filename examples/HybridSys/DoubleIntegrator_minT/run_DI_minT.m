% Double integrator with 2 modes - LQR problem
% xdot = [ x2 ] + [ 0 ] * u
%        [ 0  ]   [ 1 ]
% 
% u(t) \in [-1, 1]
% X_1 = [0.5, 2] x [-1, 1]
% X_2 = [-1, 0.5] x [-1, 1]
% G(1->2) = {0.5} x [-1, 1]\[-1e-3,1e-3]
% XT_1 = {}
% XT_2 = {(0,0)}
% h = 1
% H = 0
% 
% Trajectory starts at (0.8,0.3) in mode 1, and goes to (0,0) in minimum time
% 

clear;
T = 5;          % time horizon
d = 6;          % degree of relaxation
nmodes = 2;     % number of modes

% Define variables
t = msspoly( 't', 1 );
x = cell( nmodes, 1 );
u = cell( nmodes, 1 );
f = cell( nmodes, 1 );
g = cell( nmodes, 1 );
x0 = cell( nmodes, 1 );
hX = cell( nmodes, 1 );
hU = cell( nmodes, 1 );
hXT = cell( nmodes, 1 );
sX = cell( nmodes, nmodes );
R = cell( nmodes, nmodes );
h = cell( nmodes, 1 );
H = cell( nmodes, 1 );

x0{1} = [ 0.8; 0.3 ];

% Dynamics
x{1} = msspoly( 'x', 2 );
u{1} = msspoly( 'u', 1 );
f{1} = T * [ x{1}(2); 0 ];
g{1} = T * [ 0; 1 ];

x{2} = x{1};
u{2} = u{1};
f{2} = f{1};
g{2} = g{1};

% Domains
% Mode 1 
y = x{1};
hX{1} = [ ...
    2 - y(1);
    y(1) - 0.5;
    1 - y(2)^2 ];
hU{1} = 1 - u{1}^2;
sX{1,2} = [ ...
    y(1) - 0.5;
    0.5 - y(1);
    1 - y(2)^2;
    y(2)^2 - 1e-6 ];
R{1,2} = y;
h{1} = 1;
H{1} = 0;

% Mode 2
y = x{2};
hX{2} = [ ...
    0.5 - y(1);
    y(1) + 1;
    1 - y(2)^2 ];
hU{2} = 1 - u{2}^2;
hXT{2} = -y.^2;
h{2} = 1;
H{2} = 0;

% Options
options.freeFinalTime = 1;
options.withInputs = 1;
options.svd_eps = 1e4;

% Solve
[out] = HybridOCPDualSolver(t,x,u,f,g,hX,hU,sX,R,x0,hXT,h,H,d,options);

pval = T * out.pval;

%% Plot
xs0 = x0{1};
figure;
hold on;
box on;

% Domain
plot([0.5, 0.5], [-1, 1], 'k');

% Integrate forward trajectory
controller = [ out.u{1}; out.u{2} ];
ode_options = odeset('Events', @EventFcn);
[ tval, xval ] = ode45(@(tt,xx) T * Hybrid_DIEq( tt, xx, controller, [t;x{1}] ), ...
                       [0:0.001:1], [xs0], ode_options);
h_traj = plot(xval(:,1), xval(:,2),'LineWidth',2);

plot(x0{1}(1),x0{1}(2),'Marker','o','MarkerEdgeColor',[0 0.4470 0.7410]);
plot(0,0,'Marker','x','MarkerEdgeColor',[0 0.4470 0.7410]);
xlim([-1,2]);
ylim([-1,1]);
set(gca,'XTick',[-1,2]);
set(gca,'YTick',[-1,1]);
set(gca, 'FontSize', 20);
xlabel('$x_1$','Interpreter','LaTex','FontSize',30);
ylabel('$x_2$','Interpreter','LaTex','FontSize',30);
box on;

% Control
figure;
hold on;
uval = zeros( size(tval) );
for i = 1 : length(tval)
    if (xval(i,1) >= 0.5)
        uval(i) = double(subs(controller(1), [t;x{1}], [tval(i);xval(i,1:2)']));
    else
        uval(i) = double(subs(controller(2), [t;x{1}], [tval(i);xval(i,1:2)']));
    end
    uval(i) = min(1, max(-1,uval(i)));
end
plot(tval*T, uval,'Linewidth',2);

xlim([0,T]);


integrand = xval(:,1).^2 + xval(:,2).^2 + 20 * uval(:).^2;
cost = T * sum( integrand(1:end-1) .* diff(tval) );

disp(['Computation time = ' num2str(out.time)]);
disp(['LMI ' int2str(d) ' lower bound = ' num2str(pval)]);
disp(['Cost from simulation = ' num2str(cost)]);
