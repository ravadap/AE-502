function dydt = Delaunay_EOMS(~, y, w)
    % unpack delaunay vars
    L = y(1);

    % EOMs
    L_dot = 0;
    G_dot = 0;
    H_dot = 0;
    l_dot = 1/L^3;
    g_dot = 0;
    h_dot = w;

    dydt = [L_dot G_dot H_dot l_dot g_dot h_dot]';
    
end