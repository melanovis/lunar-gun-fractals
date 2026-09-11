format compact
clear
clc
clf reset

% ----------

TA_start_norm = 0.25;

earth_SOI_approx = 1.5e9;

G = 6.6743e-11;
error_tolerance = 1e-13;

lunar_radius = 1.7374e6;
earth_radius = 6.371e6; 

lunar_obliquity = 6.68;
earth_obliquity = 23.44;

lunar_mass = 7.35e22; %kg
earth_mass = 5.972e24;
mu = G*(earth_mass + lunar_mass);

pair_eccentricity = 0.0549;
approx_semimajor = 385692.5e3;
approx_apoapse = approx_semimajor*(1+pair_eccentricity);

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

%% setting parameters

record_video = true;

frames = 30*60;

eject_v = 2.35e3;

launch_azimuth = 90;
launch_elevation = 0;

scene_time = 3600*24*270;

launch_azimuth = rem(launch_azimuth,360);


return_struct = run_transit(lunar_start_statematrix, earth_start_statematrix, scene_time, 0, 0, 0, 0, 2.4e3, t_start, true);
earth_moon_statematrix = return_struct.earth_moon_statematrix;

scene_lunar_statematrix = [earth_moon_statematrix(:,1),earth_moon_statematrix(:,8:10)];
scene_earth_statematrix = [earth_moon_statematrix(:,1),earth_moon_statematrix(:,2:4)];

HSV_colour_bounds=[
-180,180    %hue - azimuth
1,earth_SOI_approx    %saturation - radius
0,90    %value - elevation
];

% time_range = [0,13.5*3600*24];
% yaw_bounds = [30,75];

freeze_time = 13.5*3600*24;

if record_video
    v = VideoWriter("colour_demo", 'MPEG-4');
    v.FrameRate = 30;
    v.Quality = 90;
    open(v);
end

[earth_state, ~] = interp_statematrix_via_time(scene_earth_statematrix,freeze_time);
[lunar_state, ~] = interp_statematrix_via_time(scene_lunar_statematrix,freeze_time);


start_r = [417957120  -407498816   -10075849];
[az,el,r] = cart2sph(start_r(1),start_r(2),start_r(3));
start_sph = [rad2deg(az),rad2deg(el),r];

keyframes = [
0, start_sph
0.05, -90, 0, 2e8
0.3, -90-360, 0, 2e8
0.35, -90-360+90, 0, 0
0.5, -90-360+90, 0, earth_SOI_approx
0.55, -90-360+90, 0, 2e8
0.6, -90-360+180, 0, 2e8
0.9, -90-360+180, 180, 3e8
1, start_sph(1)-360, start_sph(2), start_sph(3)
];

keyframe_time = keyframes(:,1);
keytime_duplicate = [keyframes(:,1); keyframes(:,1)-1e-10; keyframes(:,1)-2e-10];
[~,ind_sort] = sort(keytime_duplicate);
keytime_duplicate = keytime_duplicate(ind_sort);
keytime_duplicate = interp1([min(keytime_duplicate),max(keytime_duplicate)],[0,1],keytime_duplicate);

keyframe_duplicate = zeros([height(keytime_duplicate),width(keyframes)]);
keyframe_duplicate(:,1) = keytime_duplicate.';
for n=2:width(keyframes)
    keyframe_duplicate(:,n) = [keyframes(:,n); keyframes(:,n)-1e-10; keyframes(:,n)-2e-10];
    keyframe_duplicate(:,n) = keyframe_duplicate(ind_sort,n);
end

keyframe_interp = zeros([frames,width(keyframes)]);
keyframe_r = keyframe_interp;
keyframe_interp(:,1) = linspace(0,1,frames);
for n=2:width(keyframes)
    keyframe_interp(:,n) = interp1(keyframe_duplicate(:,1),keyframe_duplicate(:,n),keyframe_interp(:,1),"pchip");
end
keyframe_interp(:,1) = round(interp1([0,1],[1,frames],keyframe_interp(:,1)));

keyframe_interp(:,2) = rem(keyframe_interp(:,2)+360,360);
keyframe_interp(:,3) = rem(keyframe_interp(:,3)+360,360);

keyframe_r(:,1) = keyframe_interp(:,1);
for n=1:height(keyframe_interp)
    [key_x,key_y,key_z] = sph2cart(deg2rad(keyframe_interp(n,2)),deg2rad(keyframe_interp(n,3)),keyframe_interp(n,4));
    keyframe_r(n,2:4) = [key_x,key_y,key_z];
end


plot_horisontal_scale = 503141184;
plot_vertical_scale = 100000000;

for ind_f = 1:frames

    pixel_state = keyframe_r(ind_f,:);
    dot_RGB = get_hsv(HSV_colour_bounds, pixel_state);
    %dot_RGB = [1,1,1];

    %dot_size = 30;
    dot_size = 100;

    clf reset
    scatter(nan,nan,"k")
    hold on
    grid on
    axis tight equal

    extra(1) = plot3(earth_init_statematrix(:,2),earth_init_statematrix(:,3),earth_init_statematrix(:,4),"w");
    extra(2) = plot3(lunar_init_statematrix(:,2),lunar_init_statematrix(:,3),lunar_init_statematrix(:,4),"w");

    lunar_tick_series = linspace(0,earthluna_orbital_period,36);
    lunar_tick_series(end) = [];
    for n=1:numel(lunar_tick_series)
        [lunar_orbit_tickstate, ~] = interp_statematrix_via_time(lunar_init_statematrix,lunar_tick_series(n));
        lunar_ticks(n) = plot3([lunar_orbit_tickstate(2),lunar_orbit_tickstate(2)],[lunar_orbit_tickstate(3),lunar_orbit_tickstate(3)],[0,lunar_orbit_tickstate(4)],"--",color=[repelem(0.5,3)],LineWidth=0.03);
    end

    [surface_map_x,surface_map_y,surface_map_z] = ellipsoid(lunar_state(2),lunar_state(3),lunar_state(4),lunar_radius,lunar_radius,lunar_radius,64);
    luna_obj = surf(surface_map_x,surface_map_y,surface_map_z,edgealpha=0.1,FaceColor=[repelem(0.7,3)]);
    luna_marker = text(lunar_state(2),lunar_state(3),lunar_state(4),"$\space$ Luna",Interpreter="latex",color=[1,1,1],fontsize=14,VerticalAlignment="bottom");

    [surface_map_x,surface_map_y,surface_map_z] = ellipsoid(earth_state(2),earth_state(3),earth_state(4),earth_radius,earth_radius,earth_radius,64);
    earth_obj = surf(surface_map_x,surface_map_y,surface_map_z,edgealpha=0.1,FaceColor=[0,0,1]);
    earth_marker = text(earth_state(2),earth_state(3),earth_state(4),"$\space$ Earth",Interpreter="latex",color=[1,1,1],fontsize=14,VerticalAlignment="bottom");

    scatter3(pixel_state(2),pixel_state(3),pixel_state(4),dot_size ,"filled",MarkerFaceColor=dot_RGB,MarkerEdgeColor=[0,0,0]);
    plot3([pixel_state(2),pixel_state(2)],[pixel_state(3),pixel_state(3)],[0,pixel_state(4)],"--",color=[repelem(0.5,3)],LineWidth=0.01);
    plot3([0,pixel_state(2)],[0,pixel_state(3)],[0,pixel_state(4)],"--",color=[repelem(0.5,3)],LineWidth=0.01);
    plot3([0,pixel_state(2)],[0,pixel_state(3)],[0,0],"--",color=[repelem(0.5,3)],LineWidth=0.01);

    RGB_str = "$\quad$" + uint8(dot_RGB(1)*255) + ", " + uint8(dot_RGB(2)*255) + ", " + uint8(dot_RGB(3)*255);
    text(pixel_state(2),pixel_state(3),pixel_state(4), RGB_str,Interpreter="latex",Color=[1,1,1],FontSize=16)

    plot_horisontal_tmp = max([ plot_horisontal_scale , max(abs(pixel_state(2:3))) ]);
    plot_vertical_tmp = max([ plot_vertical_scale , abs(pixel_state(4)) ]);

    xlim([-1,1].*plot_horisontal_tmp)
    ylim([-1,1].*plot_horisontal_tmp)
    zlim([-1,1].*plot_vertical_tmp)

    fprintf("- rendering %d of %d.\n",ind_f, frames)

    view([75,25])

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

    ax.GridColor = [1 1 1];
    ax.MinorGridColor = [1 1 1];
    ax.GridAlpha = 0.2;

    legend_str = " " + newline + " " + newline+"time = " + round(freeze_time/(3600*24),2) + " days" + newline + "launch v = " + round(eject_v/1e3,4) + " km/s" + ", " + "launch lunar TA ($\nu)$ = "+round(rem(TA_start_norm*360,360),2) + "$^{\circ}$" + newline + ...
    "launch azimuth ($\theta$) = " + round(launch_azimuth,3) +"$^{\circ}$" + ", " + "launch elevation ($\epsilon$) = " + round(launch_elevation,3) +"$^{\circ}$" +" $\quad$ $\enspace$" + newline + " $\quad$";

    legend_obj = legend(legend_str,textcolor="w", FontSize=15, Interpreter="latex",location="northeast",edgecolor='k');
    %legend boxoff
    legend_obj.Color = [0,0,0]; 
    legend_obj.Units = 'normalized';
    legend_obj.Position = [0.815, 0.9, 0, 0]; 

    set(findall(gcf,'-property','FontSize'), 'FontName', 'Times')
    drawnow()
    hold off

    if record_video
        frame = [];
        frame_new = uint8([]);
        frame = getframe(gcf);
        for n=1:3
            frame_new(:,:,n) = uint8(imresize(squeeze(frame.cdata(:,:,n)), [nan, 1920*2 ]));
        end
        writeVideo(v,frame_new)
    end


end

if record_video
    close(v);
    sound(sin(2*pi*400*(0:1/14400:0.15)), 14400);
end



function sph_coords_rgb = get_hsv(HSV_colour_bounds, pixel_state)
    
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
    sph_coords_rgb = hsv2rgb(sph_coords_hsv);
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