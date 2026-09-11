format compact
clear
clc
clf reset

% ----------

TA_start_norm = 0.25;

earth_SOI_approx = 1.5e9;

G = 6.6743e-11;

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

[timerange_init, statematrix_init] = ode113(@(timerange, statematrix) earth_moon_ODE(timerange, statematrix, earth_mass, lunar_mass), timespan, earthlunar_state_intial ,odeset('Reltol',1e-13));
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
lunar_start_statematrix = interp_statematrix(lunar_init_statematrix,TA_start,8);
t_start = lunar_start_statematrix(1);
earth_start_statematrix = interp_statematrix(earth_init_statematrix,t_start,1);

%% setting parameters

record_video = true;
animation_type = 2; %1 for r, 2 for v

frames = 30*35;

eject_v = 2.45e3;

launch_azimuth = 0;
launch_elevation = 0;

scene_time = 3600*24*270;

res = ceil([720,1440]./0.75)

long_range = linspace(-180,180,res(2));
lat_range = linspace(-90,90,res(1));

lat_range_length = numel(lat_range);
long_range_length = numel(long_range);
launch_azimuth = rem(launch_azimuth,360);
pixel_structmat(1:lat_range_length, 1:lat_range_length) = struct();

if ~gcp().Connected %start up the cores
    delete(gcp('nocreate'));
    parpool('local',8);
end

parfor n=1:long_range_length
    for m=1:lat_range_length
        launch_long = long_range(n);
        launch_lat = lat_range(m);

        return_struct = run_transit_new(lunar_start_statematrix, earth_start_statematrix, scene_time, launch_long, launch_lat, launch_azimuth, launch_elevation, eject_v, t_start, 600);
        pixel_statematrix = int32(return_struct.transit_statematrix);

        switch animation_type
            case 1
                pixel_statematrix = pixel_statematrix(:,1:4);
            otherwise 
                pixel_statematrix = [pixel_statematrix(:,1),pixel_statematrix(:,5:7)];
                v_list = pixel_statematrix(:,2:4);
                v_mag = vecnorm(single(v_list), 2, 2);
                v_max = max(v_mag);
                pixel_structmat(m,n).v_norm_max = v_max;
        end

        pixel_structmat(m,n).statematrix = pixel_statematrix;
        pixel_structmat(m,n).end_time = pixel_statematrix(end,1);
    end
    fprintf("- %3.2f.\n",n/numel(long_range))
end

% save("tmp_save.mat","pixel_structmat")
% load("tmp_save.mat")

max_time = 0;
for n=1:lat_range_length
    for m=lat_range_length
        max_time = max([pixel_structmat(m,n).end_time,max_time]);
    end
end
max_time = single(max_time);

if animation_type ~= 1
    max_v = 0;
    for n=2:lat_range_length-1
        for m=lat_range_length
            max_v = max([pixel_structmat(m,n).v_norm_max,max_v]);
        end
    end
    max_v = 0.9*max_v;
    v_bounds = [0,max_v];
end

HSV_colour_bounds=[
-180,180    %hue - azimuth
1,earth_SOI_approx    %saturation - radius
0,90    %value - elevation
];

time_range = linspace(0,scene_time,frames);
% makima_factor = 3;
% interp_lower = sort(repelem(0,makima_factor) + linspace(0,1,makima_factor));
% time_range = interp1([interp_lower,max_time],[repelem(0,makima_factor),max_time],linspace(0,max_time,frames),"makima");

cmap = interp1([0,0.2,0.4,0.6,0.8,1], [[repelem(0.1,3)]; [0.259 0.039 0.408]; [0.584 0.149 0.404]; [0.867 0.318 0.227]; [0.98 0.647 0.039]; [0.98 1 0.643]], linspace(0, 1, 1e3));

if record_video
    v = VideoWriter("render_plot", 'MPEG-4');
    v.FrameRate = 30;
    v.Quality = 95;
    open(v);
end

if record_video
    w = VideoWriter("render_raw", 'MPEG-4');
    w.FrameRate = 30;
    w.Quality = 95;
    open(w);
end


for ind_t = 1:numel(time_range)

    frame_colour = uint8(zeros([size(pixel_structmat),3]));
    
    t = time_range(ind_t);
    parfor n=1:lat_range_length
        for m=1:long_range_length
            
            if t > pixel_structmat(n,m).end_time
                frame_colour(n,m,:) = [0,0,0];
                frame_greyscale(n,m) = nan;
                frame_alphamask(n,m) = 0;
            else
                
                statematrix_spec = single(pixel_structmat(n,m).statematrix);
                pixel_state = interp_statematrix(statematrix_spec,t,1);
                frame_alphamask(n,m) = 1;

                switch animation_type
                    case 1
    
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
                        
                        sph_coords_rgb = uint8(hsv2rgb(sph_coords_hsv).*255);
            
                        frame_colour(n,m,:) = sph_coords_rgb;
                    otherwise
                        frame_greyscale(n,m) = norm(pixel_state(2:4));
                end
    
    
            end
    
        end
    end
    

    if animation_type ~= 1
        frame_colour = colour_frame(frame_greyscale,cmap,v_bounds);
    end
    
    clf reset

    scatter(nan,nan,"k")
    hold on
    grid on
    axis tight equal
    fractal_map = image(frame_colour ,'xdata',[min(long_range),max(long_range)],'ydata',[min(lat_range),max(lat_range)], AlphaData = frame_alphamask);
    
    if animation_type ~= 1
        colormap(cmap)
        h = colorbar;
        set(get(h,'label'),'string','$\left| \mathrm{v} \right|$ (km/s)', Interpreter="latex", FontSize=20);
        clim(v_bounds./1e3)
        h.Color = [1,1,1];
    end
    
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
    
    legend_str = " " + newline + " " + newline+"time = " + round(time_range(ind_t)/(3600*24),2) + " days" + newline + "launch v = " + round(eject_v/1e3,4) + " km/s" + ", " + "launch lunar TA ($\nu)$ = "+round(rem(TA_start_norm*360,360),2) + "$^{\circ}$" + newline + ...
    "launch azimuth ($\theta$) = " + round(launch_azimuth,3) +"$^{\circ}$" + ", " + "launch elevation ($\epsilon$) = " + round(launch_elevation,3) +"$^{\circ}$" +" $\quad$ $\enspace$" + newline + " $\quad$";
    
    legend_obj = legend(legend_str,textcolor="w", FontSize=15, Interpreter="latex",location="northeast",edgecolor='k');
    legend_obj.Color = [0,0,0]; 
    legend_obj.Units = 'normalized';
    legend_obj.Position = [0.815, 0.95, 0, 0]; 
    
    xlabel("longitude", Interpreter="latex", FontSize=20)
    ylabel("latitude", Interpreter="latex", FontSize=20)
    
    set(findall(gcf,'-property','FontSize'), 'FontName', 'Times')

    drawnow()

    if record_video
        frame = [];
        frame_new = uint8([]);
        frame = getframe(gcf);
        for n=1:3
            frame_new(:,:,n) = uint8(imresize(squeeze(frame.cdata(:,:,n)), [nan, 1920*2 ]));
        end
        writeVideo(v,frame_new)

        writeVideo(w,flipud(frame_colour))
    end


    fprintf("rendered frame %i of %i.\n",ind_t,numel(time_range));
end

if record_video
    close(v);
    close(w);
    sound(sin(2*pi*400*(0:1/14400:0.3)), 14400);
end




function frame_colour = colour_frame(frame_greyscale, cmap, bounds)
    
    frame_colour = uint8(zeros([size(frame_greyscale),3]));
    frame_linear = reshape(frame_greyscale,1,[]);

    nan_mask = logical(isnan(frame_linear));
    
    frame_linear(nan_mask) = min(bounds);

    for n=1:3
        colour_inds = round(interp1(bounds, [1,height(cmap)], frame_linear));
        colour_inds = max(colour_inds,1);
        colour_linear = cmap(colour_inds,n);
        colour_linear(nan_mask) = 0;
        frame_spec = reshape(colour_linear,size(frame_greyscale));
        frame_colour(:,:,n) = uint8(frame_spec.*255);
    end

end


function state_out = interp_statematrix(statematrix, t, col_time)
    state_out = zeros(1,width(statematrix));
    if height(statematrix) == 1
        state_out = statematrix;
        return
    else
        state_out = zeros([1,width(statematrix)]);
        for n=1:width(statematrix)
            state_out(n) = interp1( statematrix(:,col_time),statematrix(:,n), t, "linear","extrap"); 
        end
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