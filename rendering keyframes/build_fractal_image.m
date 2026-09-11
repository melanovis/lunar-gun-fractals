function build_fractal_image(pass_matrix)

    %1. launch TA
    %2. eject v
    %3. azimuth
    %4. elevation
    %5. center x
    %6. center y
    %7. width x/2
    %8. 1440 resolution div factor
    %9. scene time

    TA_start_norm = pass_matrix(1);
    
    eject_v = pass_matrix(2);
    
    launch_azimuth = pass_matrix(3);
    launch_elevation = pass_matrix(4);
    
    scene_time = pass_matrix(9);
    
    res = ceil([720,1440]./pass_matrix(8));
    dims_x = [pass_matrix(5)-pass_matrix(7),pass_matrix(5)+pass_matrix(7)];
    dims_y = [pass_matrix(6)-pass_matrix(7)/2,pass_matrix(6)+pass_matrix(7)/2];

    long_range = linspace(dims_x(1),dims_x(2),res(2));
    lat_range = linspace(dims_y(1),dims_y(2),res(1));

    launch_azimuth = rem(launch_azimuth,360);
    
    earth_SOI_approx = 1.5e9;
    
    G = 6.6743e-11;
    error_tolerance = 1e-12;
    
    lunar_radius = 1.7374e6;
    earth_radius = 6.371e6; 
    
    lunar_obliquity = 6.68;
    earth_obliquity = 23.44;
    
    lunar_mass = 7.35e22; %kg
    earth_mass = 5.972e24;
    mu = G*(earth_mass + lunar_mass);
    
    approx_semimajor = 385692.5e3;
    
    pair_eccentricity = 0.0549;
    
    lunar_inc = 5.145;
    lunar_right_ascension = 0;
    lunar_argument_peri = 0;
    
    earth_inc = lunar_inc;
    earth_right_ascension = 0;
    earth_argument_peri = 0;
    
    %% defining earth-luna orbit
    
    %finding velocity magnitudes
    r_relative = approx_semimajor*(1-pair_eccentricity);
    v_relative = sqrt(mu*(1 + pair_eccentricity)/(approx_semimajor*(1 - pair_eccentricity)));
    
    r_earth_mag = (lunar_mass / (earth_mass + lunar_mass))*r_relative;
    r_luna_mag = (earth_mass / (lunar_mass + earth_mass))*r_relative;
    v_earth_mag = (lunar_mass / (earth_mass + lunar_mass))*v_relative;
    v_luna_mag = -(earth_mass / (earth_mass + lunar_mass))*v_relative;
    
    lunar_T_matrix = [
    cosd(lunar_right_ascension), -sind(lunar_right_ascension), 0
    sind(lunar_right_ascension), cosd(lunar_right_ascension), 0
    0, 0, 1
    ]*[
    1, 0, 0
    0, cosd(lunar_inc), -sind(lunar_inc)
    0, sind(lunar_inc), cosd(lunar_inc)
    ]*[
    cosd(lunar_argument_peri), -sind(lunar_argument_peri), 0
    sind(lunar_argument_peri), cosd(lunar_argument_peri), 0
    0, 0, 1
    ];
    
    lunar_r_initial = lunar_T_matrix*[0; r_luna_mag; 0];
    lunar_v_initial = lunar_T_matrix*[v_luna_mag; 0; 0];
    
    earth_T_matrix = [
    cosd(earth_right_ascension), -sind(earth_right_ascension), 0
    sind(earth_right_ascension), cosd(earth_right_ascension), 0
    0, 0, 1
    ]*[
    1, 0, 0
    0, cosd(earth_inc), -sind(earth_inc)
    0, sind(earth_inc), cosd(earth_inc)
    ]*[
    cosd(earth_argument_peri), -sind(earth_argument_peri), 0
    sind(earth_argument_peri), cosd(earth_argument_peri), 0
    0, 0, 1
    ];
    
    earth_r_initial = earth_T_matrix*[0; r_earth_mag; 0];
    earth_v_initial = earth_T_matrix*[v_earth_mag; 0; 0];
    
    earth_r_initial = -earth_r_initial;
    
    earthlunar_state_intial = [earth_r_initial.', earth_v_initial.', lunar_r_initial.', lunar_v_initial.'];
    
    earthluna_orbital_period = sqrt( (approx_semimajor^3)*(4*pi^2)/mu);
    
    timespan = linspace(0,earthluna_orbital_period*2,7e3);
    
    [timerange_init, statematrix_init] = ode113(@(timerange, statematrix) earth_moon_ODE(timerange, statematrix, earth_mass, lunar_mass), timespan, earthlunar_state_intial ,odeset('Reltol',error_tolerance));
    earth_init_statematrix = [timerange_init, statematrix_init(:,1:6)];
    lunar_init_statematrix = [timerange_init, statematrix_init(:,7:12)];
    
    lunar_peri_state = lunar_init_statematrix(1,:);
    earth_apo_state = earth_init_statematrix(1,:);
    
    lunar_TA_list = get_planet_TAlist(lunar_peri_state(2:4),lunar_peri_state(5:7),lunar_init_statematrix);
    earth_TA_list = get_planet_TAlist(earth_apo_state(2:4),earth_apo_state(5:7),earth_init_statematrix);
    earth_TA_list = rem(earth_TA_list+180,360); %TA function is designed to start at periapse, accounting for earth's apo start
    
    earth_TA_list = accumulate_increasing_TAlist(earth_TA_list);
    lunar_TA_list = accumulate_increasing_TAlist(lunar_TA_list);
    
    earth_init_statematrix = [earth_init_statematrix, earth_TA_list.'];
    lunar_init_statematrix = [lunar_init_statematrix, lunar_TA_list.'];
    
    
    TA_start = rem(TA_start_norm*360,360);
    lunar_start_statematrix = interp_statematrix_via_TA(lunar_init_statematrix,TA_start);
    t_start = lunar_start_statematrix(1);
    earth_start_statematrix = interp_statematrix_via_time(earth_init_statematrix,t_start);
    

    pixel_quantity = numel(long_range)*numel(lat_range);
    
    if ~gcp().Connected
        delete(gcp('nocreate'));
        parpool('local',8);
    end
    
    tic
    lat_range_length = numel(lat_range);
    parfor n=1:numel(long_range)
        for m=1:lat_range_length
            launch_long = long_range(n);
            launch_lat = lat_range(m);
            return_struct = run_transit(lunar_start_statematrix, earth_start_statematrix, scene_time, launch_long, launch_lat, launch_azimuth, launch_elevation, eject_v, t_start);
            endstate_struct(m,n) = return_struct;
            check_matrix(m,n) = logical(all(return_struct.check_matrix));
            timing_matrix(m,n) = return_struct.transit_time;
        end
        fprintf("-")
    end
    disp(toc)
    
    fprintf("\n")
    
    HSV_colour_bounds=[
    -180,180    %hue - azimuth
    1,earth_SOI_approx    %saturation - radius
    0,90    %value - elevation
    ];
    

    frame_coloured = uint8(zeros([size(endstate_struct),3]));
    
    for n=1:height(endstate_struct)
        for m=1:width(endstate_struct)
            pixel_endstate = endstate_struct(n,m).final_state;
            pixel_lifespan = endstate_struct(n,m).transit_time;
    
            pixel_state = pixel_endstate;

            black_flag = false;
            if scene_time > pixel_lifespan
                black_flag = true;
            end
   
            %spherical coord system conversion
            r_sph = [
            atan2d(pixel_state(3),pixel_state(2)) %azimuth
            norm(pixel_state(2:4)) %r
            atan2d(pixel_state(4),norm(pixel_state(2:3)) ) %elevation
            ];
        
            r_sph(3) = abs(r_sph(3) );
            
            sph_coords_hsv = [
            interp1([HSV_colour_bounds(1,:)],[0,1],r_sph(1),"linear")
            interp1([log10(HSV_colour_bounds(2,:))],[0,1],log10(r_sph(2)),"linear")
            interp1([HSV_colour_bounds(3,:)],[1,0],r_sph(3),"linear")
            ].';
    
            %sph_coords_hsv(3) = 1;
            
            sph_coords_rgb = hsv2rgb(sph_coords_hsv).*255;
    
            if black_flag
                sph_coords_rgb = [0,0,0];
            end
    
            frame_coloured(n,m,:) = uint8(sph_coords_rgb);    
        end
    end

    scatter(nan,nan,"k")
    hold on
    grid on
    axis tight equal
    xlim([min(long_range),max(long_range)])
    ylim([min(lat_range),max(lat_range)])
    plotmap = frame_coloured;
    wm = image(plotmap ,'xdata',[min(long_range),max(long_range)],'ydata',[min(lat_range),max(lat_range)]);

    ax = gca;
    ax.FontSize = 20;

    set(gcf, 'Color', [0,0,0])
    set(ax, 'Color', [0,0,0])

    ax.XColor = 'w';
    ax.YColor = 'w';
    ax.ZColor = 'w';
    ax.Title.Color  = 'w';
    ax.XLabel.Color = 'w';
    ax.YLabel.Color = 'w';

    xlabel("longitude", Interpreter="latex", FontSize=20)
    ylabel("latitude", Interpreter="latex", FontSize=20)

    legend_str = " " + newline + " " + newline+"time = " + round(scene_time/(3600*24),2) + " days" + newline + "launch v = " + round(eject_v/1e3,4) + " km/s" + ", " + "launch lunar TA ($\nu)$ = "+round(rem(TA_start_norm*360,360),2) + "$^{\circ}$" + newline + ...
    "launch azimuth ($\theta$) = " + round(launch_azimuth,3) +"$^{\circ}$" + ", " + "launch elevation ($\epsilon$) = " + round(launch_elevation,3) +"$^{\circ}$" +" $\quad$ $\enspace$" + newline + " $\quad$";

    legend_obj = legend(legend_str,textcolor="w", FontSize=15, Interpreter="latex",location="northeast",edgecolor='k');
    legend_obj.Color = [0,0,0]; 
    legend_obj.Units = 'normalized';
    legend_obj.Position = [0.815, 0.87, 0, 0]; 
    
    set(findall(gcf,'-property','FontSize'), 'FontName', 'Times')
    drawnow()

end


function state_out = interp_statematrix_via_TA(statematrix,TA_target)
    row_time = 8;
    [~,inds_between] = mink(abs(statematrix(:,row_time)-TA_target),2);
    inds_between = sort(inds_between);
    for n=1:width(statematrix)
        state_out(n) = interp1( [statematrix(inds_between(1),row_time),statematrix(inds_between(2),row_time)], [statematrix(inds_between(1),n),statematrix(inds_between(2),n)], TA_target);
    end
end

function [state_out, inds_between] = interp_statematrix_via_time(statematrix,t_target)
    row_time = 1;
    [~,inds_between] = mink(abs(statematrix(:,row_time)-t_target),2);
    inds_between = sort(inds_between);
    for n=1:width(statematrix)
        state_out(n) = interp1( [statematrix(inds_between(1),row_time),statematrix(inds_between(2),row_time)], [statematrix(inds_between(1),n),statematrix(inds_between(2),n)], t_target,"linear","extrap");
    end
end



function TA_list_out = accumulate_increasing_TAlist(TA_list_in)
    a = 0;
    TA_prev = TA_list_in(1);
    for n=1:numel(TA_list_in)
        TA_curr = TA_list_in(n);
        if abs(TA_curr - TA_prev) > 180 %assume this means a rollback from 360 to 0
            a = a+360;
        end
        TA_list_out(n) = TA_curr + a;
        TA_prev = TA_curr;
    end
end

function TA_list = get_planet_TAlist(a,b,statematrix)
    a = a/norm(a);
    b = b/norm(b);

    plane_normal = cross(a,b);

    if plane_normal(3) < 0
        plane_normal = -plane_normal;
    end

    for n=1:height(statematrix)
        c = statematrix(n,2:4);
        TA_list(n) = rem( atan2d( dot(plane_normal,cross(a,c)) , dot(a,c) ) + 360, 360);
    end
end


function state_out = earth_moon_ODE(t, state_in, earth_mass, lunar_mass)

    G = 6.6743e-11;

    r_earth = state_in(1:3).';
    v_earth = state_in(4:6).';
    r_luna = state_in(7:9).';
    v_luna = state_in(10:12).';

    r_earthluna = r_earth - r_luna;
    r_unit_earthluna = norm(r_earthluna);

    acc_earth =  -G * lunar_mass * r_earthluna / r_unit_earthluna^3;
    acc_luna = G * earth_mass * r_earthluna / r_unit_earthluna^3;

    %state in derived for the next timestep
    state_out = [v_earth, acc_earth, v_luna, acc_luna].';
end

function state_out = earth_twobody_ODE(t,state,primary_mass)
    G = 6.67e-11;
    r = state(1:3); 
    v = state(4:6);
    r_mag = norm(r);
    r_unit = r/r_mag;
    force_gravity = r_unit*(-G*primary_mass/(r_mag^2));
    acc = force_gravity;
    state_out = [v; acc];
end