function xi = kepler_universal_anomaly(mu, dt, r0, vr0, alpha)
    
    % Initial estimate of u_anomaly
    xi = sqrt(mu)*abs(alpha)*dt;

    % Set tolerance and initialize ratio
    tol = 1e-8;
    ratio = 1;

    % Newton's Method
    while abs(ratio) > tol
        % Calculate f(universal anomaly) = f(xi)
        zi = alpha*xi^2;
        f_xi = ((r0*vr0)/(sqrt(mu)))*xi^2*C(zi) + (1-alpha*r0)*xi^3*S(zi) + r0*xi - sqrt(mu)*dt;

        % Calculate f'(xi)
        df_xi = ((r0*vr0)/(sqrt(mu)))*xi*(1-alpha*xi^2*S(zi)) + (1-alpha*r0)*xi^2*C(zi) + r0;

        % Calculate the ratio f(xi)/f'(xi)
        ratio = f_xi/df_xi;

        % Update value of xi
        xi = xi - ratio;
    end

end