function [r,v] = r_v(mu, dt, r0, v0)
    
    % magnitude of r0 and v0
    r0_mag = norm(r0);
    v0_mag = norm(v0);

    % radial component of velocity
    vr0 = dot(r0, v0)/r0_mag;

    % alpha
    alpha = 2/r0_mag - v0_mag^2/mu;

    % find universal anomaly
    x = kepler_universal_anomaly(mu, dt, r0_mag, vr0, alpha);

    % f and g
    [f,g] = f_and_g(mu, alpha, r0_mag, dt, x);

    % r and r_mag
    r = f*r0 + g*v0;
    r_mag = sqrt(dot(r,r));

    % fdot and gdot
    [fdot, gdot] = fdot_and_gdot(mu, alpha, r0_mag, r_mag, x);

    % v
    v = fdot*r0 + gdot*v0;
end