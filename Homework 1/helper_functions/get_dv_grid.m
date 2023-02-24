function [dv_1I_r,dv_1I_f] = get_dv_grid(dt_arrival, dt_departure, mu_sun, r_traj_arrival, r_traj_departure, v_traj_arrival, v_traj_departure, traj_type, TU)
    
    dv_1I_r = NaN(length(dt_arrival), length(dt_departure));
    dv_1I_f = dv_1I_r;

    for i=1:length(dt_arrival) % arrival
        for j=1:length(dt_departure) % departure
            % [i,j]
            TOF = days(dt_arrival(i)-dt_departure(j))/(TU/86400);
            if TOF > 0
                [v1,v2] = lambert_solver(mu_sun, r_traj_departure(j,:), r_traj_arrival(i,:), TOF, traj_type);
                if ~isnan(v1(1)) && ~isnan(v2(1))
                    dv_1I_r(i,j) = norm(v1 - v_traj_departure(j,:)) + norm(v_traj_arrival(i,:) - v2);
                    dv_1I_f(i,j) = norm(v1 - v_traj_departure(j,:));
                end
            end
        end
    end

end