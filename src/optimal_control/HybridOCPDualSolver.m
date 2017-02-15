function [out] = HybridOCPDualSolver(t,x,u,f,g,hX,hU,sX,R,x0,hXT,h,H,d,options)
% Hybrid optimal control problem - dual solver
% ------------------------------------------------------------------------
% t     -- time indeterminate,      1-by-1 free msspoly
% x     -- state indeterminate,     I-by-1 cell of (n_i)-by-1 free msspoly
% u     -- control indeterminate,   I-by-1 cell of (m_i)-by-1 free msspoly
% f     -- dynamics, part 1,        I-by-1 cell of (n_i)-by-1 msspoly in x{i}
% g     -- dynamics, part 2,        I-by-1 cell of (n_i)-by-(m_i) msspoly in x{i}
% hX    -- domain (X_i),            I-by-1 cell of (~)-by-1 msspoly in x{i}
% hU    -- set U_i,                 I-by-1 cell of (~)-by-1 msspoly in u{i}
% sX    -- guard(i,j),              I-by-I cell of (~)-by-1 msspoly in x{i}
% R     -- reset map(i,j),          I-by-I cell of (n_j)-by-1 msspoly in x{i}
% x0    -- initial point,           I-by-1 cell of (n_i)-by-1 reals
% hXT   -- target set,              I-by-1 cell of (~)-by-1 msspoly in x{i}
% h     -- running cost,            I-by-1 cell of 1-by-1 msspoly in (t,x{i},u{i})
% H     -- terminal cost,           I-by-1 cell of 1-by-1 msspoly in x{i}
% d     -- degree of relaxation,    positive even number (scalar)
% options  -- struct
% ------------------------------------------------------------------------
% Solves the following hybrid optimal control problem:
% 
% inf  { \int_0^1 h(t,x,u) dt } + H( x(1) )
% s.t. \dot{x_i} = f_i + g_i * u
%      x_j(t+) = R_ij( x_i(t-) ) when { x_i(t-) \in S_ij }
%      x(0) = x0
%      x(t) \in X
%      x(T) \in XT
%      u(t) \in U
% 
% where:
% X_i  = { x | each hX{i}(x) >= 0 }          domain of mode i
% U_i  = { u | each hU{i}(u) >= 0 }          range of control in mode i
% S_ij = { x | each sX{i}(x) >= 0 } <= X_i   guard in mode i
% XT_i = { x | each hXT{i}(x)>= 0 } <= X_i   target set in mode i
% 
% We solve the following *weak formulation* via relaxation of degree 'd':
% 
% sup  v_i(0,x_0)
% s.t. LFi(v_i) + h_i >= 0
%      v_i(T,x) <= H_i
%      v_i(t,x) <= v_j( t, R_ij (x) )
% ------------------------------------------------------------------------
% 'options' is a struct that contains:
%     .freeFinalTime:   1 = free final time, 0 = fixed final time (default: 0)
%     .MinimumTime:     1 = minimum time, 0 = o.w. (default: 0)
%     .withInputs:      1 = perform control synthesis, 0 = o.w. (default: 0)
%     .solver_options:  options that will be passed to SDP solver (default: [])
%     .svd_eps:         svd threshold (default: 1e3)
% ------------------------------------------------------------------------
% Attention: T = 1 is fixed number. To solve OCP for different time
% horizons, try scaling the dynamics and cost functions.
% 

%% Sanity check
if mod(d,2) ~= 0
    warning('d is not even. Using d+1 instead.');
    d = d+1;
end
nmodes = length(x);

max_m = 0;
for i = 1 : nmodes
    m = length(u{i});
    if m > max_m
        max_m = m;
    end
    if (length(f{i}) ~= length(x{i})) || (size(g{i},1) ~= length(x{i}))
        error('Inconsistent matrix size.');
    end
    if size(g{i},2) ~= m
        error('Inconsistent matrix size.');
    end
end

svd_eps = 1e3;
if nargin < 15, options = struct(); end
if isfield(options, 'svd_eps'), svd_eps = options.svd_eps; end
if ~isfield(options, 'freeFinalTime'), options.freeFinalTime = 0; end
if ~isfield(options, 'withInputs'), options.withInputs = 0; end

T = 1;
hT = t*(T-t);

if isempty(R)
    R = cell(nmodes,nmodes);
    for i = 1 : nmodes
    for j = 1 : nmodes
        if ~isempty(sX{i,j}), R{i,j} = x{i}; end
    end
    end
end

%% Setup spotless program
% define the program
prog = spotsosprog;
prog = prog.withIndeterminate( t );

% constraints to keep track of
mu_idx = zeros(nmodes, 1);

% Create program variables in each mode
for i = 1 : nmodes
    
    prog = prog.withIndeterminate( x{i} );
    prog = prog.withIndeterminate( u{i} );
    
    % create v(i)
    vmonom{ i } = monomials( [ t; x{ i } ], 0:d );
    [ prog, v{ i }, ~ ] = prog.newFreePoly( vmonom{ i } );
    
    % create the variables that will be used later
    vT{ i } = subs( v{ i }, t, T );
    dvdt{ i } = diff( v{ i }, t );
    dvdx{ i } = diff( v{ i }, x{ i } );
    Lfv{ i } = dvdt{ i } + dvdx{ i } * f{ i };
    Lgv{ i } = dvdx{ i } * g{ i };
    Lv{ i } = Lfv{ i } + Lgv{ i } * u{ i };
end

% creating the constraints and cost function
obj = 0;
for i = 1 : nmodes
    % Lv_i + h_i >= 0                   Dual: mu
    prog = sosOnK( prog, Lv{ i } + h{ i }, ...
                   [ t; x{ i }; u{ i } ], [ hT; hX{ i }; hU{ i } ], d);
    mu_idx(i) = size( prog.sosExpr, 1 );
    
    % v(T,x) <= H_i(x)                  Dual: muT
    if ~isempty( hXT{ i } )
        if options.freeFinalTime
            prog = sosOnK( prog, H{ i } - v{ i }, [ t; x{ i } ], [ hT; hXT{ i } ], d );
        else
            prog = sosOnK( prog, H{ i } - vT{ i }, x{ i }, hXT{ i }, d );
        end
    end
    
    % v_i(t,x) <= v_j (t,R(x))          Dual: muS
    for j = 1 : nmodes
        if ( ~isempty( sX{ i, j } ) ) % if its empty there isn't a guard between these
            vj_helper = subs( v{ j }, x{ j }, R{ i, j } ); 
            prog = sosOnK( prog, vj_helper - v{ i }, [ t; x{ i } ] , [ hT; sX{ i, j } ], d );
        end
    end
    
    % Objective function
    if ~isempty( x0{ i } )
        obj = obj + subs( v{ i }, [ t; x{i} ], [ 0; x0{i} ] );
    end
end

% set options
spot_options = spot_sdp_default_options();
spot_options.verbose = 1;

if isfield(options, 'solver_options')
    spot_options.solver_options = options.solver_options;
end

%% Solve
tic;
[sol, y, dual_basis, ~] = prog.minimize( -obj, @spot_mosek, spot_options );
out.time = toc;

out.pval = double(sol.eval(obj));
out.sol = sol;

%% Control Synthesis
if ~options.withInputs
    out.u = [];
    return;
end

uout = cell( nmodes, max_m );
u_real_basis = cell( nmodes, 1 );
for i = 1 : nmodes
    urb = monomials( [ t; x{ i } ], 0 : d/2 );
    u_real_basis{ i } = urb;
    % moments of mu
    mu_basis = dual_basis{ mu_idx(i) };
    y_mu = sol.dualEval( y{ mu_idx(i) } );
    
    % moment matrix of mu
    mypoly = mss_s2v( urb * urb' );
    coeff_mu = DecompMatch( mypoly, mu_basis );
    moment_mat = double( mss_v2s( coeff_mu * y_mu ) );
    
    [U,S,V] = svd(moment_mat);
    iS1 = S;

    startExp = 1e-10;
    while norm(pinv(iS1)) / norm(S) > svd_eps
        iS1(iS1 < startExp) = 0;
        startExp = startExp * 10;
    end

    fprintf('norm of moment matrix %0.2f\n', norm(S));
    fprintf('norm of inverse moment matrix %0.2f\n', norm(pinv(iS1)));

    iMyt = V*pinv(iS1)*U';
    
    % yp and yn
    for j = 1 : length( u{ i } )
        fprintf('Processing mode %1.0f, input #%1.0f ...\n', i, j );
        mypoly = u{ i }( j ) * urb;
        coeff = DecompMatch( mypoly, mu_basis );
        y_u = sol.dualEval( coeff * y_mu );
        
        u_coeff = iMyt * y_u;
        uout{ i, j } = u_coeff' * urb;
        
    end
end

out.u = uout;
