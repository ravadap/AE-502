function F_z = F(mu, r1_mag, r2_mag, dt, A, z)
    F_z = (y(r1_mag, r2_mag, A, z)/C(z))^(3/2)*S(z)+A*sqrt(y(r1_mag, r2_mag, A, z))-sqrt(mu)*dt;
end