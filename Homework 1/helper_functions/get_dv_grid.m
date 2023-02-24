function dv = get_dv_grid(dt_arrival, dt_departure, mu_sun, r_traj_arrival, r_traj_departure, v_traj_arrival, v_traj_departure, traj_type)
    
    dv = zeros(length(dt_arrival), length(dt_departure));

    for i=1:length(dt_arrival) % arrival
        for j=1:length(dt_departure) % departure
            [i,j]
            TOF = dt_arrival(i)-dt_departure(j);
            if TOF > 0
                % Check minimum TOF
                % tp = min_TOF(r_traj_e1(j,:), r_traj_1I(i,:), mu_sun, traj_type);
                % if TOF >= tp
                    [v1,v2] = lambert_solver(mu_sun, r_traj_departure(j,:), r_traj_arrival(i,:), TOF, traj_type)
                    if ~isnan(v1(1)) && ~isnan(v2(1))
                        dv(i,j) = norm(v1 - v_traj_departure(j,:)) + norm(v_traj_arrival(i,:) - v2);
                    end
                % end
            end
        end
    end

end