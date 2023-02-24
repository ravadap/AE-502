function [v1,v2]=lambert_solver(mu, r1, r2, dt, trajectory)
    
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

    % Calculate A
    A = sin(df)*sqrt((r1_mag*r2_mag)/(1-cos(df)));

    % Set tolerance and ratio
    tol = 1e-10;
    ratio = 1;
    z = 0;

    val = 0;
    while val <= 0
        z = z + 1;
        val = F(mu, r1_mag, r2_mag, dt, A, z);
    end

  
    iter = 0;
    while abs(ratio) > tol
        iter = iter + 1;
        F_z = F(mu, r1_mag, r2_mag, dt, A, z);
        dF_z = (1/(2*sqrt(y(r1_mag, r2_mag, A, z)*C(z)^5)))*((2*C(z)*dS(z)-3*dC(z)*S(z))*y(r1_mag, r2_mag, A, z)^2+(A*C(z)^(5/2)+3*C(z)*S(z)*y(r1_mag, r2_mag, A, z))*dy(A,z));
        ratio = F_z/dF_z;
        z = z-ratio;

        if iter > 1000
            v1 = nan;
            v2 = v1;
            return
        end
    end

    z = real(z);
    
    % Calculate y
    y_z = y(r1_mag, r2_mag, A, z); 

    % Calculate f, g, and gdot functions
    f = 1-y_z/r1_mag;
    g = A*sqrt(y_z/mu);
    gdot = 1-y_z/r2_mag;

    % Calculate v1 and v2
    v1 = (1/g)*(r2-f*r1);
    v2 = (1/g)*(gdot*r2-r1);

%     if iter > 100
%         v1 = [0,0,0];
%         v2 = v1;
%     end

end