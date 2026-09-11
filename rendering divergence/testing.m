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
lunar_start_statematrix = interp_statematrix(lunar_init_statematrix,TA_start,8);
t_start = lunar_start_statematrix(1);
earth_start_statematrix = interp_statematrix(earth_init_statematrix,t_start,1);

%% setting parameters


eject_v = 2.35e3;

launch_azimuth = 90;
launch_elevation = 0;

scene_time = 3600*24*270;

launch_azimuth = rem(launch_azimuth,360);

res = ceil([720,1440]./3)
%res = ceil([960,1920]./1.5)

long_range = linspace(-180,180,res(2));
lat_range = linspace(-90,90,res(1));

pixel_quantity = numel(long_range)*numel(lat_range);

if ~gcp().Connected %start up the cores
    delete(gcp('nocreate'));
    parpool('local',8);
end


animation_type = 3; %1 - long peturb only, 2 - lat peturb only, 3 - both peturb
frames = 100;
run_sim = true;

peturb_factor = 5e-7;

lat_range_length = numel(lat_range);
long_range_length = numel(long_range);

if run_sim
    pixel_structmat(1:lat_range_length, 1:lat_range_length) = struct();
    parfor n=1:long_range_length
        for m=1:lat_range_length
    
            switch animation_type
                case 1
                    launch_long = long_range(n)-peturb_factor;
                    launch_lat = lat_range(m);
                case 2
                    launch_long = long_range(n);
                    launch_lat = lat_range(m)-peturb_factor;
                otherwise 
                    launch_long = long_range(n)-peturb_factor;
                    launch_lat = lat_range(m)-peturb_factor;
            end
    
            return_struct_peturbed = run_transit(lunar_start_statematrix, earth_start_statematrix, scene_time, launch_long, launch_lat, launch_azimuth, launch_elevation, eject_v, t_start);
            statematrix_peturbed_back = return_struct_peturbed.transit_statematrix;
    
            switch animation_type
                case 1
                    launch_long = long_range(n)+peturb_factor;
                    launch_lat = lat_range(m);
                case 2
                    launch_long = long_range(n);
                    launch_lat = lat_range(m)+peturb_factor;
                otherwise 
                    launch_long = long_range(n)+peturb_factor;
                    launch_lat = lat_range(m)+peturb_factor;
            end
    
            return_struct_peturbed = run_transit(lunar_start_statematrix, earth_start_statematrix, scene_time, launch_long, launch_lat, launch_azimuth, launch_elevation, eject_v, t_start);
            statematrix_peturbed_forward = return_struct_peturbed.transit_statematrix;
    

            state_dr_norm = zeros([height(statematrix_peturbed_back),1]);
            state_dv_norm = zeros([height(statematrix_peturbed_back),1]);

            for p = 1:height(statematrix_peturbed_back)
                state_peturbed_back = statematrix_peturbed_back(p,1);
                state_peturbed_forward = interp_statematrix(statematrix_peturbed_forward, statematrix_peturbed_back(p,1) ,1);

                state_delta = single( (state_peturbed_forward - state_peturbed_back) / (2*peturb_factor) );
                %state_delta = single( (state_peturbed_forward - state_peturbed_back) );

                state_dr_norm(p,1) = norm(state_delta(2:4));
                state_dv_norm(p,1) = norm(state_delta(5:7));
            end

            state_dr_norm = round(state_dr_norm);
            state_dv_norm = round(state_dv_norm);

            pixel_structmat(m,n).time = statematrix_peturbed_back(:,1)';
            pixel_structmat(m,n).max_time = statematrix_peturbed_back(end,1);
            % pixel_structmat(m,n).state_dr_norm = state_dr_norm;
            % pixel_structmat(m,n).state_dv_norm = state_dv_norm;
            pixel_structmat(m,n).state_dr_norm = movmean(state_dr_norm, ceil(height(statematrix_peturbed_back)/50));
            pixel_structmat(m,n).state_dv_norm = movmean(state_dv_norm, ceil(height(statematrix_peturbed_back)/50));
    
            % plot(state_dr_norm)
            % set(gca,"yscale","log")
            % drawnow()
            % pause

        end
        fprintf("- %3.2f.\n",n/numel(long_range))
    end
    save("tmp_struct_save.mat","pixel_structmat");
end


load("tmp_struct_save.mat");

%bounds
max_dr = 0;
max_dv = 0;
max_time = 0;
for n=1:long_range_length
    for m=lat_range_length
        max_dr = max([max(pixel_structmat(m,n).state_dr_norm),max_dr]);
        max_dv = max([max(pixel_structmat(m,n).state_dv_norm),max_dv]);
        max_time = max([max(pixel_structmat(m,n).max_time),max_time]);
    end
end

cmap = interp1([0,0.2,0.4,0.6,0.8,1], [[repelem(0.1,3)]; [0.259 0.039 0.408]; [0.584 0.149 0.404]; [0.867 0.318 0.227]; [0.98 0.647 0.039]; [0.98 1 0.643]], linspace(0, 1, 1e3));

timerange = linspace(0,max_time,frames);

ind_frame = 100;

frame_grey_dr = zeros([lat_range_length,long_range_length]);
t = timerange(ind_frame);

parfor n=1:height(pixel_structmat)
    for m=1:long_range_length

        if t > pixel_structmat(n,m).max_time
            %out of bounds or collided
            frame_grey_dr(n,m) = nan;
            frame_grey_dv(n,m) = nan;
        else
            dr_spec = interp1(pixel_structmat(n,m).time, pixel_structmat(n,m).state_dr_norm, t, "makima");
            frame_grey_dr(n,m) = dr_spec;

            dv_spec = interp1(pixel_structmat(n,m).time, pixel_structmat(n,m).state_dv_norm, t, "makima");
            frame_grey_dv(n,m) = dv_spec;
        end
    end
end
frame_grey_dr([1,end],:) = nan;
frame_grey_dv([1,end],:) = nan;

colormap(jet)

subplot(2,1,1)
hold on
grid on
axis tight equal
imagesc(log10(frame_grey_dr))
%set(gca,"colorscale","log")
%clim([0,max_dr])

subplot(2,1,2)
hold on
grid on
axis tight equal
imagesc(log10(frame_grey_dv))
%set(gca,"colorscale","log")
%clim([0,max_dv])



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