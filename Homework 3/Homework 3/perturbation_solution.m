function [kepler_state, rv_state] = perturbation_solution(mu, init_vals, tol, perturb)
    % unpack initial values
    a = init_vals(1);
    e = init_vals(2);
    i = init_vals(3);
    l0 = init_vals(4);
    g0 = init_vals(5);
    h0 = init_vals(6);

    % calculate other initial values
    L0 = sqrt(a);
    G0 = L0 * sqrt(1-e^2);
    H0 = G0 * cosd(i);

    state = [L0 G0 H0 l0 g0 h0];

    % use ODE45 to integrate from t0 to tf
    t0 = 0;
    tf = 100; % 100 TU
    tspan = linspace(t0, tf, 10000);
    options = odeset('RelTol',tol, 'AbsTol',tol);
    
    [t,y] = ode89(@Delaunay_EOMS, tspan, state, options, perturb);
    
    L = y(:,1);
    G = y(:,2);
    H = y(:,3);
    l = y(:,4); % mean anomaly
    g = y(:,5); % argument of perigee
    h = y(:,6); % argument of perigee + RAAN
    
    % Convert to Keplerian 
    a = L.^2;
    e = sqrt(1-(G./L).^2);
    i = acos(H./G); % rad
    M = l;
    w = g;
    OM = h-g;
    
    f = zeros(size(M)); % in rad
    for j=1:length(M)
        f(j) = get_f_from_M(M(j), e(j), tol);
    end

    kepler_state = [a, e, i, w, OM, f];
    
    rt = zeros(length(tspan), 3);
    vt = zeros(length(tspan), 3);
    
    for j=1:length(tspan)
        [rt(j,:), vt(j,:)] = elm2rv_PR(a(j), e(j), i(j), w(j), OM(j), f(j), mu);
    end

    rv_state = [rt, vt];
end