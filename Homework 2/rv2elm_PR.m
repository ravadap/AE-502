% Convert functions from state vectors to keplerian orbital elements

function [a,e_mag,i,w,OM,f] = rv2elm(mu, r,v)
    % Calculate the magnitude of the position and velocity
    r_mag = norm(r);
    v_mag = norm(v);

    % Radial velocity
    vr_mag = dot(r,v)/r_mag;

    % Angular momentum
    h = cross(r,v);
    h_mag = norm(h);

    % inclination
    i = acos(h(3)/h_mag);

    % RAAN Right ascencion of ascending node
    k_hat = [0 0 1];
    N = cross(k_hat, h);
    N_mag = norm(N);

    if N(2) >= 0
        OM = acos(N(1)/N_mag);
    else
        OM = 2*pi - acos(N(1)/N_mag);
    end

    % Eccentricity
    e = (1/mu)*((v_mag^2-mu/r_mag)*r-r_mag*vr_mag*v);
    e_mag = norm(e);
   
    % Argument of Perigee
    if e(3) >= 0
        w = acos(dot(N, e)/(N_mag*e_mag));
    else
        w = 2*pi - acos(dot(N, e)/(N_mag*e_mag));
    end

    % True anomaly
    if vr_mag >= 0 
        f = acos(dot(e,r)/(e_mag*r_mag));
    else
        f = 2*pi - acos(dot(e,r)/(e_mag*r_mag));
    end

    % Semi-major axis
    rp = (h_mag^2/mu)*(1/(1+e_mag*cosd(0)));
    ra = (h_mag^2/mu)*(1/(1+e_mag*cosd(180)));
    a = 0.5*(rp+ra);

end