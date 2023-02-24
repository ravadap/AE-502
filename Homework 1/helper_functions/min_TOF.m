function tp = min_TOF(r1, r2, mu, trajectory)
    % Calculate r1 and r2
    r1_mag = norm(r1);
    r2_mag = norm(r2);

    % Choose prograde or retrograde trajectory and calculate df
    r_cross = cross(r1,r2);

    switch trajectory

        % Prograde trajectory (0 < i < 90)
        case 'prograde'
            if r_cross(3) >= 0
                df = acos(dot(r1,r2)/(r1_mag*r2_mag));
            else
                df = 2*pi - acos(dot(r1,r2)/(r1_mag*r2_mag));
            end

        % Retrograde trajectory (90 < i < 180)
        case 'retrograde'
            if r_cross(3) < 0
                df = acos(dot(r1,r2)/(r1_mag*r2_mag));
            else
                df = 2*pi - acos(dot(r1,r2)/(r1_mag*r2_mag));
            end
    end

    c = sqrt(r1_mag^2 + r2_mag^2 - 2*r1_mag*r2_mag*cos(df));
    s = (r1_mag+r2_mag+c)/2;
    
    tp = sqrt(2)/(3*sqrt(mu))*(s^(3/2)-(sign(sin(df))*(s-c)^(3/2)));
end