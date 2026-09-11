function return_struct = run_transit(lunar_start_statematrix, earth_start_statematrix, scene_time, launch_long, launch_lat, launch_azimuth, launch_elevation, eject_v, t_start)

    warning("off","all")

    earth_mass = 5.972e24;
    lunar_mass = 7.35e22; 

    collision_tolerance = 10e3; 

    lunar_inc = 5.145;
    lunar_obliquity = 6.68;
    lunar_radius = 1.7374e6;
    earth_radius = 6.371e6; 
    earth_SOI_approx = 1.5e9;

    error_tolerance = 5e-6;

    perp_lunar_plane = cross(lunar_start_statematrix(2:4),[0,0,1]);
    perp_lunar_plane = perp_lunar_plane./norm(perp_lunar_plane);
    
    pole_unshifted = [0,0,1];
    pole_theta = lunar_inc - lunar_obliquity;
    pole_corrected = pole_unshifted*cosd(pole_theta) + cross(perp_lunar_plane,pole_unshifted)*sind(pole_theta);
    
    null_point = pole_corrected*cosd(90) + cross(perp_lunar_plane,pole_corrected)*sind(90); %this is 0,0 long/lat
    equator_90 = cross(pole_corrected,null_point);
    
    pole_norm = pole_corrected./norm(pole_corrected);
    null_norm = null_point./norm(null_point);
    equator_90_norm = equator_90./norm(equator_90);
    
    if abs(launch_lat) == 90
        launch_lat = sign(launch_lat)*(90-1e-4);
    end
    
    launch_r = (cosd(launch_lat)*cosd(launch_long))*null_norm + (cosd(launch_lat)*sind(launch_long))*equator_90_norm + sind(launch_lat)*pole_norm;
    
    launch_r_norm = (launch_r./norm(launch_r));
    launch_r_lunar = launch_r_norm * (lunar_radius+10);
    
    launch_tang_plane = pole_norm - dot(pole_norm,launch_r_norm)*launch_r_norm;
    launch_ref_north = launch_tang_plane./norm(launch_tang_plane);
    launch_ref_east = cross(launch_r_norm, launch_ref_north);
    
    tang_vector = cosd(-launch_azimuth)*launch_ref_north + sind(-launch_azimuth)*launch_ref_east;
    eject_vector = cosd(launch_elevation)*tang_vector + sind(launch_elevation)*launch_r_norm;
    
    eject_vector_norm = eject_vector/norm(eject_vector);
    
    eject_v_geo = eject_vector_norm*eject_v + lunar_start_statematrix(5:7);
    eject_r_geo = lunar_start_statematrix(2:4) + launch_r_lunar;
    
    total_transit_statematrix = [];
    
    
    scene_timerange = [];
    scene_statematrix = [];
    rocket_start_statematrix = [t_start, eject_r_geo, eject_v_geo];
    ER3BP_state_inital = [rocket_start_statematrix(2:7), earth_start_statematrix(2:7), lunar_start_statematrix(2:7)];
    ER3BP_timespan = [0,scene_time];
    ode_settings = odeset('RelTol', error_tolerance,'Events', @boundary_check);
    [scene_timerange, scene_statematrix] = ode113(@(scene_timerange, scene_statematrix) threebody_ODE(scene_timerange, scene_statematrix, earth_mass, lunar_mass), ER3BP_timespan, ER3BP_state_inital , ode_settings);
    scene_timerange = scene_timerange + t_start;
    scene_eject_statematrix = [scene_timerange, scene_statematrix]; %statematrix before first burn
    
    total_transit_statematrix = scene_eject_statematrix;
    
    lunar_approach = sqrt(sum((total_transit_statematrix(:,2:4) - total_transit_statematrix(:,14:16)).^2, 2));
    ind_lunarcoltol = 1;
    while ind_lunarcoltol < numel(lunar_approach)
        if lunar_approach(ind_lunarcoltol) > collision_tolerance + lunar_radius
            break
        else
            ind_lunarcoltol = ind_lunarcoltol+1;
        end
    end
    lunar_approach_secondary = lunar_approach;
    lunar_groundingcheck = lunar_approach_secondary(1:ind_lunarcoltol-1);
    lunar_approach_secondary(1:ind_lunarcoltol-1) = inf;
    
    [~,ind_min_l1] = min(lunar_approach);
    lunar_r_closest = norm(total_transit_statematrix(ind_min_l1,2:4) - total_transit_statematrix(ind_min_l1,14:16));
    
    [~,ind_min_l2] = min(lunar_approach_secondary);
    lunar_r_closest_secondary = norm(total_transit_statematrix(ind_min_l2,2:4) - total_transit_statematrix(ind_min_l2,14:16)); 
    
    earth_approach = sqrt(sum((total_transit_statematrix(:,2:4) - total_transit_statematrix(:,8:10)).^2, 2));
    [~,ind_min_e] = min(earth_approach);
    earth_r_closest = norm(total_transit_statematrix(ind_min_e,2:4) - total_transit_statematrix(ind_min_e,8:10));
    
    collision_flag = false;
    if lunar_r_closest < lunar_radius || earth_r_closest < (earth_radius + collision_tolerance) || lunar_r_closest_secondary < (lunar_radius + collision_tolerance) || any(lunar_groundingcheck < lunar_radius)
        collision_flag = true;
    end
    
    worst_collision_intersect = 0;
    if collision_flag
        worst_collision_intersect = min([lunar_r_closest - lunar_radius, lunar_r_closest_secondary - (lunar_radius + collision_tolerance), earth_r_closest - (earth_radius + collision_tolerance) ]);
    end
    
    transit_time = total_transit_statematrix(end,1) - total_transit_statematrix(1,1);
    
    height_original = height(total_transit_statematrix);
    
    if collision_flag
        if lunar_r_closest < lunar_radius
            ind_collision = ind_min_l1;
        elseif lunar_r_closest_secondary < (lunar_radius + collision_tolerance)
            ind_collision = ind_min_l2;
        else 
            ind_collision = ind_min_e;
        end
    
        total_transit_statematrix(ind_collision+1:end,:)=[];
    end
    
    deep_space_flag = false;
    if max(earth_approach) > earth_SOI_approx && ~collision_flag
        deep_space_flag = true;
        [~,SOI_exit_crossing] = mink(abs(earth_approach-earth_SOI_approx),10);
        ind_SOI_exit = min(SOI_exit_crossing);
        total_transit_statematrix(ind_SOI_exit+1:end,:) = [];
    end

    check_matrix = [
    ~collision_flag
    ~deep_space_flag
    ~any(any(isnan(total_transit_statematrix)))
    ];

    %total_transit_statematrix(:,5:end) = [];
    total_transit_statematrix = round(total_transit_statematrix,3);

    total_transit_statematrix(:,1) = total_transit_statematrix(:,1) - total_transit_statematrix(1,1);

    return_struct = struct();
    return_struct = setfield(return_struct,"final_state",single(total_transit_statematrix(end,:)));
    return_struct = setfield(return_struct,"check_matrix",logical(check_matrix));
    return_struct = setfield(return_struct,"transit_time",single(total_transit_statematrix(end,1)));
end

function [distance, terminate, return_state] = boundary_check(t,state)

    r_barycenter = norm(state(1:3));
    r_earth = norm(state(1:3) - state(7:9));
    r_luna = norm(state(1:3) - state(13:15));

    distance = [
    r_barycenter - (1.5e9+1e7)
    r_earth - (6.371e6+5e3+1)
    r_luna - (1.7374e6)
    ];

    terminate = [
    1
    1
    1
    ];

    return_state=0;
end