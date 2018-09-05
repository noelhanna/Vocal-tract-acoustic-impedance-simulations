
function [ output ] = ...
    Vocal_Tract_General_Model_2017_SVT(tube, glottis, SG, PARAMS, Flex)
% calculates the impedance etc. of a series of cylinders
% use with .mat file that specifies tube.length, tube.radius
% and tube.alpha_multiplier
% can choose to specify glottis, subglottal tract and PARAMS
% note glottis is always a single tube but tube and SG can be concatenated
% otherwise default values are used


%% Define default parameters (if none specified)

if nargin < 5
    disp('Default NH non-rigid tract properties used')
    % compliance of tissue
    Flex.C_tissue = 10.9*10^-8;
%     Flex.C_tissue = 2.9*10^-8; % male, 4.9*10^-8; %   female value
    % inertance of tissue
    Flex.L_tissue = 1400 ;
%     Flex.L_tissue = 1400 ; % male & female value
    % damping coeffiecient
    % female - resistance due to low frequency (real part of impedance)
    
    Flex.R_tissue = 1*10^5;
    
    % Flex.R_tissue = ((1.8*10^5) + (10*10^5) )/2;
%     Flex.R_tissue = 0.86*10^5; % male % 10*10^5; % female - use bA0 only
    % Flex.R_tissue = ((2.7*10^5) + (8.6*10^5) )/2; %5.5*10^5; %9.4*10^5;
    % average male value%1.6*10^5; % 4.2*10^4;
    % resistance due to low frequency (real part of impedance)
    
    if nargin < 4
        disp('default speed of sound parameters used')
%         PARAMS.Temp = 37;
        PARAMS.Temp = 20.2;
        PARAMS.humidity = 0.4;
        % used for the visco-thermal loss calculations
        PARAMS.DeltaT = PARAMS.Temp-26.85;
        
        PARAMS.resolution = 0.01; % frequency resolution in Hz
        PARAMS.freq = 0.000001:PARAMS.resolution:4000;
        PARAMS.w = 2.*pi.*PARAMS.freq; % angular frequency
        
        Parameters.temperature = PARAMS.Temp;
        Parameters.humidity = PARAMS.humidity;
%         PARAMS.c = 359; % 37 degrees in moist air
        PARAMS.c = 340; % standard
%         PARAMS.c = calculatespeedofsound(Parameters);
        PARAMS.rho = 1.14; % air density
        
        if nargin < 3
            disp('default subglottal tract used')
            SG.length(1,1) = 0.20;
            SG.radius(1,1) = 0.01;
            SG.alpha_multiplier(1,1) = 1;
%             SG.length(1,1) = 0.0000001;
%             SG.length(2,1) = 1.2*sum(tube.length);
%             SG.radius(1,1) = 0.9*nanmean(tube.radius);
%             SG.radius(2,1) = 0.9*nanmean(tube.radius);
%             SG.alpha_multiplier(1,1) = 1;
%             SG.alpha_multiplier(2,1) = 1;
            
            if nargin < 2
                disp('default glottis used')
                glottis.length = 0.01;
                % NB if radius < 0.00001 then it will be ideally closed
                glottis.radius = 0.000001;
%                 glottis.radius = 0.01;
                glottis.alpha_multiplier = 1;
                
                if nargin < 1
                    disp('default single cylinder model used')
                    tube.length = 0.17;
                    tube.radius = 0.0131;
                    tube.alpha_multiplier = 5;
                end
            end
        end
    end
end

tube.glottis = glottis;
tube.SG = SG;

%% define the calibration (radiation) impedance
% seen by the source for a single microphone ACUZ measurement

sourcetube.length = 10000;
sourcetube.radius = 0.003;
sourcetube.alpha_multiplier = 1;

% characteristic impedance
sourcetube.Z0(1, :) = z0(PARAMS, sourcetube, 1);

% wavenumber with adjusted attenuation coefficent
sourcetube.K(1, :) = wavenumber(PARAMS, sourcetube, 1);

% flange load impedance
tube.Z_rad_source(1, :) = zradiation(sourcetube, 1);

%% define flange load impedance at lips
% (the lip end is the highest numbered tube segment)
lip_counter = length(tube.length);

% characteristic impedance
tube.Z0(lip_counter, :) = z0(PARAMS, tube, lip_counter);

% wavenumber with adjusted attenuation coefficent
tube.K(lip_counter, :) = wavenumber(PARAMS, tube, (lip_counter));

% flange load impedance
tube.Z_rad(lip_counter, :) = zradiation(tube, lip_counter);

tube.Z_load(lip_counter, :) = tube.Z_rad(lip_counter, :);

% input impedance with flange termination
tube.Z(lip_counter, :) ...
    = zinput(tube, tube.Z_load(lip_counter,:), lip_counter);

%%%%%%%%%%% define outgoing wave amplitude as 1 - this is temporary
tube.A(lip_counter, :) = ones(1,length(PARAMS.freq));

% calculate pressure and flow at output and input of each segment
[tube.p_out(lip_counter,:) ...
    tube.p_in(lip_counter,:) ...
    tube.u_out(lip_counter,:) ...
    tube.u_in(lip_counter,:) ...
    tube.BonA(lip_counter, :)] ...
    = calc_pressure_flow(tube, lip_counter);

% Power = |p|*|u|*cos(angle(p-u))
tube.P_out(lip_counter,:) ...
    = abs(tube.p_out(lip_counter,:))...
    .*abs(tube.u_out(lip_counter,:))...
    .*cos( angle(tube.p_out(lip_counter,:)...
    -tube.u_out(lip_counter,:) ) );

tube.P_in(lip_counter,:) ...
    = abs(tube.p_in(lip_counter,:))...
    .*abs(tube.u_in(lip_counter,:))...
    .*cos( angle(tube.p_in(lip_counter,:)...
    -tube.u_in(lip_counter,:) ) );

% use a simple function for transfer functions TF = x/y
% use TF = transfer_fn(specified_out, specified_in)


%% General recursion for all other tube segments
for k = lip_counter-1 : -1 : 1 % count backwards from inside lips
    % characteristic impedance
    tube.Z0(k,:) = z0(PARAMS, tube, k);
    
    % wavenumber with adjusted attenuation coefficent
    tube.K(k,:) = wavenumber(PARAMS, tube, k);
    
    tube.Z_load(k, :) = tube.Z(k+1, :);
    
    % load impedance is the next tube, i.e. tube.Z(k+1)
    tube.Z(k,:) = zinput(tube, tube.Z(k+1,:), k);
    
    % assume outgoing wave amplitude is 1 - temporary
    tube.A(k, :) = ones(1,length(PARAMS.freq));
    
    
    % % % % % calculate all transfer functions
tube.T_pp(lip_counter,:) = ...
    ( ( abs(tube.p_out(lip_counter,:)) )...
    ./( abs(tube.p_in(lip_counter,:)) ) );
tube.T_pu(lip_counter,:) = ...
    ( ( abs(tube.p_out(lip_counter,:)) )...
    ./( abs(tube.u_in(lip_counter,:)) ) );
tube.T_uu(lip_counter,:) = ...
    ( ( abs(tube.u_out(lip_counter,:)) )...
    ./( abs(tube.u_in(lip_counter,:)) ) );
tube.T_PP(lip_counter,:) = ...
    ( ( real(tube.P_out(lip_counter,:)) )...
    ./( real(tube.P_in(lip_counter,:)) ) );
tube.T_Pp(lip_counter,:) = ...
    ( ( real(tube.P_out(lip_counter,:)) )...
    ./( abs(tube.p_in(lip_counter,:)) ) );
tube.T_Pu(lip_counter,:) = ...
    ( ( real(tube.P_out(lip_counter,:)) )...
    ./( abs(tube.u_in(lip_counter,:)) ) );
    
    % calculate pressure and flow at output and input of segment
    [tube.p_out(k,:) ...
        tube.p_in(k,:) ...
        tube.u_out(k,:) ...
        tube.u_in(k,:) ...
        tube.BonA(k,:)] ...
        = calc_pressure_flow(tube, k);
    
    % Power = |p|*|u|*cos(theta)
    tube.P_out(k,:) ...
        = abs(tube.p_out(k,:))...
        .*abs(tube.u_out(k,:))...
        .*cos( angle(tube.p_out(k,:)...
        -tube.u_out(k,:) ) );
    
    tube.P_in(k,:) ...
        = abs(tube.p_in(k,:))...
        .*abs(tube.u_in(k,:))...
        .*cos( angle(tube.p_in(k,:)...
        -tube.u_in(k,:) ) );
    
    % %     % calculate all transfer functions, N.B. each one is
    % %     % multiplide by the previous k+1th transfer function
    tube.T_pp(k,:) = ...
        ( ( abs(tube.p_out(k,:)) )./( abs(tube.p_in(k,:)) ) )...
        .*tube.T_pp(k+1,:);
    tube.T_pu(k,:) = ...
        ( ( abs(tube.p_out(k,:)) )./( abs(tube.u_in(k,:)) ) )...
        .*tube.T_pu(k+1,:);
    tube.T_uu(k,:) = ...
        ( ( abs(tube.u_out(k,:)) )./( abs(tube.u_in(k,:)) ) )...
        .*tube.T_uu(k+1,:);
    tube.T_PP(k,:) = ...
        ( ( real(tube.P_out(k,:)) )./( real(tube.P_in(k,:)) ) )...
        .*tube.T_PP(k+1,:);
    tube.T_Pp(k,:) = ...
        ( ( real(tube.P_out(k,:)) )./( abs(tube.p_in(k,:)) ) )...
        .*tube.T_Pp(k+1,:);
    tube.T_Pu(k,:) = ...
        ( ( real(tube.P_out(k,:)) )./( abs(tube.u_in(k,:)) ) )...
        .*tube.T_Pu(k+1,:);
    
end

%% include non-rigid properties

Z_tissue_comp = zcompliance(PARAMS.w, Flex.C_tissue);
Z_tissue_inert = zinert(PARAMS.w, Flex.L_tissue);

% add the inertance and compliance in series
Ztissue = Z_tissue_inert + Z_tissue_comp;
% add the damping coefficient
Ztissue = Ztissue + Flex.R_tissue;

%%%%%%%%%%%% need to scale the amount added to the radius of the tube %%
for o = 1:lip_counter
    % add non-rigid impedance to walls in parallel to tube
    tube.Z_flex(o,:) = zparallel(tube.Z(o, :), Ztissue);
end

% calculate pressure and flow for lip (last) segment
[tube.p_out_flex(lip_counter,:) ...
    tube.p_in_flex(lip_counter,:) ...
    tube.u_out_flex(lip_counter,:) ...
    tube.u_in_flex(lip_counter,:) ...
    tube.BonA(lip_counter, :)] ...
    = calc_pressure_flow(tube, (lip_counter));

% Power = |p||u|*cos(theta)
tube.P_out_flex(lip_counter,:) ...
    = abs(tube.p_out_flex(lip_counter,:))...
    .*abs(tube.u_out_flex(lip_counter,:))...
    .*cos( angle(tube.p_out_flex(lip_counter,:)...
    -tube.u_out_flex(lip_counter,:) ) );

tube.P_in_flex(lip_counter,:) ...
    = abs(tube.p_in_flex(lip_counter,:))...
    .*abs(tube.u_in_flex(lip_counter,:))...
    .*cos( angle(tube.p_in_flex(lip_counter,:)...
    -tube.u_in_flex(lip_counter,:) ) );

% % % % % calculate all transfer functions
tube.T_pp_flex(lip_counter,:) = ...
    ( ( abs(tube.p_out_flex(lip_counter,:)) )...
    ./( abs(tube.p_in_flex(lip_counter,:)) ) );
tube.T_pu_flex(lip_counter,:) = ...
    ( ( abs(tube.p_out_flex(lip_counter,:)) )...
    ./( abs(tube.u_in_flex(lip_counter,:)) ) );
tube.T_uu_flex(lip_counter,:) = ...
    ( ( abs(tube.u_out_flex(lip_counter,:)) )...
    ./( abs(tube.u_in_flex(lip_counter,:)) ) );
tube.T_PP_flex(lip_counter,:) = ...
    ( ( real(tube.P_out_flex(lip_counter,:)) )...
    ./( real(tube.P_in_flex(lip_counter,:)) ) );
tube.T_Pp_flex(lip_counter,:) = ...
    ( ( real(tube.P_out_flex(lip_counter,:)) )...
    ./( abs(tube.p_in_flex(lip_counter,:)) ) );
tube.T_Pu_flex(lip_counter,:) = ...
    ( ( real(tube.P_out_flex(lip_counter,:)) )...
    ./( abs(tube.u_in_flex(lip_counter,:)) ) );


% recursive calculations for all other segments
for n = lip_counter-1:-1:1
    % flexible (non-rigid) pressure and flow
    [tube.p_out_flex(n,:) ...
        tube.p_in_flex(n,:) ...
        tube.u_out_flex(n,:) ...
        tube.u_in_flex(n,:) ...
        tube.BonA(n, :)] ...
        = calc_pressure_flow(tube, n);
    
    % Power = pu*cos(theta)
    tube.P_out_flex(n,:) ...
        = abs(tube.p_out_flex(n,:))...
        .*abs(tube.u_out_flex(n,:))...
        .*cos( angle(tube.p_out_flex(n,:)...
        -tube.u_out_flex(n,:) ) );
    
    tube.P_in_flex(n,:) ...
        = abs(tube.p_in_flex(n,:))...
        .*abs(tube.u_in_flex(n,:))...
        .*cos( angle(tube.p_in_flex(n,:)...
        -tube.u_in_flex(n,:) ) );
    
    % % % % % calculate all transfer functions, N.B. each one is multiplied by
    % % % the previous n+1th transfer function
    tube.T_pp_flex(n,:) = ...
        ( ( abs(tube.p_out_flex(n,:)) )...
        ./( abs(tube.p_in_flex(n,:)) ) ).*tube.T_pp_flex(n+1,:);
    tube.T_pu_flex(n,:) = ...
        ( ( abs(tube.p_out_flex(n,:)) )...
        ./( abs(tube.u_in_flex(n,:)) ) ).*tube.T_pu_flex(n+1,:);
    tube.T_uu_flex(n,:) = ...
        ( ( abs(tube.u_out_flex(n,:)) )...
        ./( abs(tube.u_in_flex(n,:)) ) ).*tube.T_uu_flex(n+1,:);
    tube.T_PP_flex(n,:) = ...
        ( ( real(tube.P_out_flex(n,:)) )...
        ./( real(tube.P_in_flex(n,:)) ) ).*tube.T_PP_flex(n+1,:);
    tube.T_Pp_flex(n,:) = ...
        ( ( real(tube.P_out_flex(n,:)) )...
        ./( abs(tube.p_in_flex(n,:)) ) ).*tube.T_Pp_flex(n+1,:);
    tube.T_Pu_flex(n,:) = ...
        ( ( real(tube.P_out_flex(n,:)) )...
        ./( abs(tube.u_in_flex(n,:)) ) ).*tube.T_Pu_flex(n+1,:);
    
end


%% use lungs as Z_load (impedance in the other direction)
SG.Z0(1,:) = z0(PARAMS, SG, 1);
SG.K(1,:) = wavenumber(PARAMS, SG, 1);

% flange load impedance
tube.Z_lung = zradiation(SG, 1);

% SG tract
tube.Z_SG(1,:) = zinput(SG, tube.Z_lung, 1);

% if more than 1 tube of SG tract is specified
if length(SG.radius) > 1
for k = 2:1:length(SG.radius) % count from lungs to glottis
    % characteristic impedance
    SG.Z0(k,:) = z0(PARAMS, SG, k);
    
    % wavenumber with adjusted attenuation coefficent
    SG.K(k,:) = wavenumber(PARAMS, SG, k);
    
    % load impedance is the next tube, i.e. tube.Z(k-1)
    tube.Z_SG(k,:) = zinput(SG, tube.Z_SG(k-1,:), k);
end
end

glottis.Z0 = z0(PARAMS, glottis, 1);
glottis.K = wavenumber(PARAMS, glottis, 1);

% -----------------------------------------
% TEMP 2017 use Flex params in parallel with SG
Z_SG_comp = zcompliance(PARAMS.w, Flex.C_tissue);
Z_SG_inert = zinert(PARAMS.w, Flex.L_tissue);

% add the inertance and compliance in series
ZSGflex = Z_SG_inert + Z_SG_comp;
% add the damping coefficient
ZSGflex = ZSGflex + Flex.R_tissue;
tube.Z_SG_flex = zparallel(tube.Z_SG(end, :), ZSGflex);
% -----------------------------------------

% glottis
tube.Z_glottis = zinput(glottis, tube.Z_SG(end,:), 1);

if glottis.radius < 0.00001;
% use ideally stopped glottis
tube.Z_in(1, :) = zstop( tube, 1);
else
% use glottal impedance loaded by SG tract
tube.Z_in(1, :) = zinput(tube, tube.Z_glottis, 1);
end

%%%%%%%%%%%%% add non-rigidity to walls in parallel to tube %%%%%%
%%%%%%% perhaps not same parameters%%
tube.Z_in_flex(1,:) = zparallel(tube.Z_in(1, :), Ztissue);


%% recursion for inward impedance
for m = 2:lip_counter
    % replace with load impedance of next tube, i.e. tube.Z(m-1)
    tube.Z_in(m,:) ...
        = zinput(tube, tube.Z_in(m-1, :), m);
    
    %%%%%%%%%%%%%% check how this is calculated - scale with radius? %%%%%
    tube.Z_in_flex(m,:) ...
        = zparallel(tube.Z_in(m, :), Ztissue);
    
end

%% Z and Transfer functions from glottis

% rigid VT
tube.out.Z = tube.Z(1,:);
tube.out.T_pp = abs(tube.p_out(lip_counter,:))./abs(tube.p_in(1,:));
tube.out.T_pu = abs(tube.p_out(lip_counter,:))./abs(tube.u_in(1,:));
tube.out.T_uu = abs(tube.u_out(lip_counter,:))./abs(tube.u_in(1,:));
tube.out.T_PP = real(tube.P_out(lip_counter,:))./real(tube.P_in(1,:));
tube.out.T_Pp = real(tube.P_out(lip_counter,:))./abs(tube.p_in(1,:));
tube.out.T_Pu = real(tube.P_out(lip_counter,:))./abs(tube.u_in(1,:));


%%%%%%% adjust the realtive amplitude of each tube section
% initial amplitude of outgoing wave = 1, i.e. 100%
tube.A(1, :) = ones(1,length(PARAMS.freq));

%%% calculate adjusted amplitude assuming that the outgoing wave initially
% had amplitude A = 1. Each segment then reduces A by BonA
for mm = 2:lip_counter
    tube.A(mm,:) = adjust_amplitude(tube.A(mm-1, :), tube.BonA(mm-1, :));
end


for k = 1:lip_counter % count backwards from inside lips
    % recalculate pressure, flow and power with the adjusted wave amplitudes
    
    % calculate pressure and flow at output and input of segment
    [tube.p_out_A(k,:), ...
        tube.p_in_A(k,:), ...
        tube.u_out_A(k,:), ...
        tube.u_in_A(k,:), ...
        tube.BonA_A(k,:)] ...
        = calc_pressure_flow(tube, k);
    
    % Power = |p|*|u|*cos(theta)
    tube.P_out_A(k,:) ...
        = abs(tube.p_out_A(k,:))...
        .*abs(tube.u_out_A(k,:))...
        .*cos( angle(tube.p_out_A(k,:)...
        -tube.u_out_A(k,:) ) );
    
    tube.P_in_A(k,:) ...
        = abs(tube.p_in_A(k,:))...
        .*abs(tube.u_in_A(k,:))...
        .*cos( angle(tube.p_in_A(k,:)...
        -tube.u_in_A(k,:) ) );
    
end


% transfer functions adjusted for wave transmission
tube.out.T_pp_A = abs(tube.p_out_A(lip_counter,:))...
    ./abs(tube.p_in_A(1,:));
tube.out.T_pu_A = abs(tube.p_out_A(lip_counter,:))...
    ./abs(tube.u_in_A(1,:));
tube.out.T_uu_A = abs(tube.u_out_A(lip_counter,:))...
    ./abs(tube.u_in_A(1,:));
tube.out.T_PP_A = real(tube.P_out_A(lip_counter,:))...
    ./real(tube.P_in_A(1,:));
tube.out.T_Pp_A = real(tube.P_out_A(lip_counter,:))...
    ./abs(tube.p_in_A(1,:));
tube.out.T_Pu_A = real(tube.P_out_A(lip_counter,:))...
    ./abs(tube.u_in_A(1,:));



% non-rigid VT
tube.out.Z_flex = tube.Z_flex(1,:);
tube.out.T_pp_flex = ...
    abs(tube.p_out_flex(lip_counter,:))./abs(tube.p_in_flex(1,:));
tube.out.T_pu_flex = ...
    abs(tube.p_out_flex(lip_counter,:))./abs(tube.u_in_flex(1,:));
tube.out.T_uu_flex = ...
    abs(tube.u_out_flex(lip_counter,:))./abs(tube.u_in_flex(1,:));
tube.out.T_PP_flex = ...
    real(tube.P_out_flex(lip_counter,:))...
    ./real(tube.P_in_flex(1,:));
tube.out.T_Pp_flex = ...
    real(tube.P_out_flex(lip_counter,:))...
    ./abs(tube.p_in_flex(1,:));
tube.out.T_Pu_flex = ...
    real(tube.P_out_flex(lip_counter,:))...
    ./abs(tube.u_in_flex(1,:));


%% Z and T from lips
tube.in.Z = tube.Z_in(lip_counter,:);

example_C = 0.01*10^-8;
Z_example_comp  = zcompliance( PARAMS.w, example_C );
tube.in.Z_with_compliance = zparallel(tube.in.Z, Z_example_comp);

tube.in.Z_flex = tube.Z_in_flex(lip_counter,:);


%ratios
tube.in.Z_ratio_1 = (zparallel(tube.Z_rad(end,:), tube.in.Z))...
    ./(tube.Z_rad(end,:));
tube.in.Z_ratio_2 = (zparallel(tube.Z_rad(end,:), tube.in.Z))...
    ./(tube.Z_rad_source);
tube.in.Z_ratio_3 = (zparallel(tube.Z_rad_source, tube.in.Z))...
    ./(tube.Z_rad_source(end,:));
tube.in.Z_ratio_4 = zparallel(tube.Z_rad_source,...
    (zparallel(tube.Z_rad(end,:), tube.in.Z)))...
    ./(tube.Z_rad(end,:));



%% Output data and calculate bandwidths
% % % formant = formant_bandwidth(PARAMS, tube.T_pp_flex(1,:));
% % tube.formant.f = formant.max.f;
% % tube.formant.z = formant.max.z;
% % tube.formant.b = formant.max.b;


impedance = impedance_bandwidth(PARAMS, tube.Z(1,:));
tube.impedance.max_f = impedance.max.f;
tube.impedance.max_z = impedance.max.z;
tube.impedance.max_b = impedance.max.b;

tube.impedance.min_f = impedance.min.f;
tube.impedance.min_z = impedance.min.z;
tube.impedance.min_b = impedance.min.b;


% calculate tube area and radius as a function of length
% length_ = [tube.SG.length; tube.glottis.length+sum(tube.SG.length);...
%     sum(tube.glottis.length)+sum(tube.SG.length)+tube.length];
% % 
% for g = 2:length(tube.length)
%     length_(g+2) = tube.length(g)+length_(g+1);
%     g
% end

length_SG(1) = tube.SG.length(1);
for g = 2:length(SG.length)
    length_SG(g) = SG.length(g)+length_SG(g-1);
end

length_glottis(1) = length_SG(end)+glottis.length(1);
for h = 2:length(glottis.length)
    length_glottis(h) = glottis.length(h)+length_glottis(h-1);
end

length_tube(1) = length_glottis(end)+tube.length(1);
for j = 2:length(tube.length)
    length_tube(j) = tube.length(j)+length_tube(j-1);
end
    
length_ = [length_SG'; length_glottis'; length_tube'];

radius_ = [tube.SG.radius; tube.glottis.radius; tube.radius];

area_ = (radius_.^2).*pi;

tube.area_fn = [length_ area_];
tube.radius_fn = [length_ radius_];

% adjust profile to glottal 0
tube.distance = tube.radius_fn(:,1)-sum(tube.SG.length)...
    -sum(tube.glottis.length);

% plot(tube.area_fn(:,1), tube.area_fn(:,2), 's-'); hold all


output.freq = PARAMS.freq;

tube.Flex = Flex;

output.tube = tube;

output.PARAMS = PARAMS;

end









%% Background functions
% The functions below are called from the main code above



function soundSpeed = calculatespeedofsound(Parameters)
% calculatespeedofsound Calculates the speed of sound
% based on Owen Cramer (1993) J. Acoust. Soc. Am. 93(5)
% p2510-2616; formula at p. 2514
% SPEED_SOUND(t,h) calculates the speed of sound given
% the temperature in Celsius t and the relative humidity h 
% (expressed as a fraction).
% The speed is calculated at one atm pressure and typical mole fraction 
% of CO_2.
if (Parameters.humidity < 0) || (Parameters.humidity > 1),...
        error('The relative humidity must be between 0 and 1.'); end
% Use p = 1 atm and C0_2 mole fraction given by Cramer Table 1 (p. 2511)
PRESSURE_ATMOSPHERIC = 101325;
X_C = 0.000314;
% Calculate mole fraction of water vapour using equation in 
% Cramer Appendix
% (p. 2515)
moleFractionOfWater = 1.00062 + 3.14e-8*PRESSURE_ATMOSPHERIC...
    + 5.6e-7 * Parameters.temperature^2;
temperatureKelvin = Parameters.temperature + 273.15;
p_sv = exp(1.2811805e-5*temperatureKelvin^2 - 1.9509874e-2...
    * temperatureKelvin + 34.04926034 - 6.3536311e3 / temperatureKelvin);
x_w = Parameters.humidity * moleFractionOfWater ...
    * p_sv/PRESSURE_ATMOSPHERIC;
soundSpeed = calculatespeedofsoundcramer(Parameters,...
    PRESSURE_ATMOSPHERIC,x_w,X_C);
end

function soundSpeed = calculatespeedofsoundcramer(Parameters,...
    PRESSURE_ATMOSPHERIC,x_w,X_C)
% SPEED_SOUND_CRAMER Calculate the speed of sound.
% Based on Owen Cramer (1993) J. Acoust. Soc. Am. 93(5) p2510-2616; 
% formula at p2514
% SPEED_SOUND_CRAMER(t,p,x_w,x_c) calculates the speed of sound given the
% temperature in Celsius t, the pressure in Pa p and the mole fraction of
% water vapour and CO_2 x_w and x_c.
if (Parameters.temperature < 0) || (Parameters.temperature > 30),...
        error('The temperature must be between 0 and 30 degrees C.'); end
if (PRESSURE_ATMOSPHERIC < 75e3) || (PRESSURE_ATMOSPHERIC > 102e3),...
        error ('The pressure must be between 75 and 102 kPa.'); end
if (x_w < 0) || (x_w > 0.06),...
        error ('The H2O mole fraction must be between 0 and 0.06.'); end
if (X_C < 0) || (X_C > 0.01),...
        error ('The CO2 mole fraction must be between 0 and 0.01.'); end
a = [
    331.5024
    0.603055
    -0.000528
    51.471935
    0.1495874
    -0.000782
    -1.82e-7
    3.73e-8
    -2.93e-10
    -85.20931
    -0.228525
    5.91e-5
    -2.835149
    -2.15e-13
    29.179762
    0.000486
    ]';
coeff = [
    1
    Parameters.temperature
    Parameters.temperature^2
    x_w
    Parameters.temperature*x_w
    Parameters.temperature^2*x_w
    PRESSURE_ATMOSPHERIC
    Parameters.temperature*PRESSURE_ATMOSPHERIC
    Parameters.temperature^2*PRESSURE_ATMOSPHERIC
    X_C
    Parameters.temperature*X_C
    Parameters.temperature^2*X_C
    x_w^2
    PRESSURE_ATMOSPHERIC^2
    X_C^2
    x_w*PRESSURE_ATMOSPHERIC*X_C
    ];
soundSpeed = a*coeff;
end


%% characteristic impedance of a cylinder of specified radius
function [ Z0 ] = z0(PARAMS, tube, counter)
rho = PARAMS.rho;
c = PARAMS.c;
radius = tube.radius(counter);

Z0 = rho*c/(pi.*(radius.^2));
end

%% impedance of an inertance with length l, area A,
function [ Zinert ] = zinert( w, L )
Zinert = 1i.*w.*L;

end

%% impedance of a known compliance C = Volume/rho*c^2
function [ Zcomp ] = zcompliance( w, C )
% Zcompliance = 1/jwC,
Zcomp = 1./(1i.*w.*C);

end


%% complex wavenumber including wall losses
function [ K ] = wavenumber(PARAMS, tube, counter)

freq = PARAMS.freq;
w = PARAMS.w;
c = PARAMS.c;
DeltaT = PARAMS.DeltaT;

radius = tube.radius(counter);
alpha_multiplier = tube.alpha_multiplier(counter);

rv = 632.8...
    .*radius...
    .*(freq.^0.5)...
    .*(1-0.0029...
    .*(DeltaT)); % approx

rt = 532.2...
    .*radius...
    .*(freq.^0.5)...
    .*(1-0.0031...
    .*(DeltaT)); % approx

% Cp./Cv ratio of specific heats (approx 1.403 at 0, 1.401 at 100)
gamma = 1.400;

v = c...
    .*( 1 - (1./(rv.*sqrt(2)))...
    -((gamma-1)./(rt.*sqrt(2))) );

% attenuation coefficient alpha
alpha = (w./c)...
    .*( (1./(rv.*sqrt(2))) + ((gamma-1)./(rt.*sqrt(2))) );

% adjust alpha to approximate real VT measurements (alpha_multiplier 5)
new_alpha = alpha.*alpha_multiplier;

% complex propogation coefficient
K = ((w./v) - (1i.*new_alpha));

end

%% flanged opening radiation impedance
function [ Z_flange ] = zradiation(tube, counter)
% RADIATION IMPEDANCE CALCULATION FROM DALMONT ET AL. 2001
% calculates the radiation impedance with a flange

K = tube.K(counter,:);
radius = tube.radius(counter);
Z0 = tube.Z0(counter);

ka = K*radius;

% % % % simple approximation from Fletcher and Rossing
% % % % real part of radiation impedance (F&R eq 8.29)
% % % Z_flange_R = Z0.*( (((ka).^2)./2) - (((ka).^4)./(2.*2.*3))...
% % %     + (((ka).^6)./(2.*2.*3.*3.*4)) );
% % % % imaginary part of radiation impedance (F&R eq 8.30)
% % % Z_flange_X = (Z0./(pi.*ka.^2)).*( (((2.*ka).^3)./3)...
% % %     - (((2.*ka).^5)./(3.*3.*5)) ...
% % %     + (((2.*ka).^7)./(3.*3.*5.*5.*7)) );
% % % 
% % % Z_flange = Z_flange_R + 1i.*Z_flange_X;


% alternative approximation from Dalmont et al. (2001)
% define the end correction for the low frequency limit
d_simple = 0.8216 * radius;

% determine the frequency-dependent end correction (15a)
d_simple = d_simple ./ (1 + (0.77 * ka).^2 ./ (1 + 0.77 * ka));

% determine the modulus of the reflection coefficient (15b)
modR = (1 + (0.323 .* ka) - (0.077 .* ka.^2))...
    ./ (1 + (0.323 .* ka) + ((1 - 0.077) .* ka.^2));

% calculate the imaginary part of the end correction
% since R = -e^(-2kjd(complex)) = -modR*e^(-2kjd(real))
% so, d(complex) = d*ln(modR))
di = log(modR)./(2 * K);
d = d_simple + 1i.*di;

% calculate the impedance (9)
Z_flange = 1i * Z0 * tan(K .* d); % Flanged opening

end

%% Input impedance of a cylinder with known Z0, Z_load, K and L
function [ Zin ] = zinput(tube, Z_load, counter)

Z0 = tube.Z0(counter);
K = tube.K(counter,:);
L = tube.length(counter);

% from Fletcher and Rossing Chapter 8 (8.23)
Zin = Z0...
    .*( ( (Z_load.*cos(K.*L)) + (1i.*Z0.*sin(K.*L)) ) ...
    ./ ( (1i.*Z_load.*sin(K.*L)) + (Z0.*cos(K.*L)) ) );
end

%% impedance of an ideally stopped end Z = -jZocotkL
function [ Zstop ] = zstop( tube, counter)

K = tube.K(counter,:);
L = tube.length(counter);
Z0 = tube.Z0(counter);

% from Fletcher and Rossing Chapter 8 (8.24)
% Zstop = -jZocotkL
Zstop = -1i.*Z0.*cot(K.*L);

end

%% Parallel addition
function [ Zparallel ] = zparallel(Z1, Z2)
% parallel addition of impedances
Zparallel = (Z1.*Z2)./(Z1 + Z2);
end


%% determine input and output pressure and flow at each segment
function [ p_out, p_in, u_out, u_in, BonA ] = ...
    calc_pressure_flow(tube, counter)
% from J Smith via Fletcher and Rossing

K = tube.K(counter,:);
L = tube.length(counter);
Z_load = tube.Z_load(counter,:);
Z0 = tube.Z0(counter);
A = tube.A(counter);

% exponential term
eTerm = exp(1i.*K.*L);

% definitions p = A/eTerm + BeTerm, u = 1/Z0 (A/eTerm - BeTerm)
% reflection coefficient is B/A
BonA = (Z_load - Z0)...
    ./((Z_load + Z0).*eTerm.*eTerm);

% pressure calculations
p_out = A.*((1./eTerm) + (BonA.*eTerm));
p_in = A.*(1 + BonA);
% flow calculations
u_out = A.*(1./Z0).*((1./eTerm) - (BonA.*eTerm));
u_in = A.*(1./Z0).*(1 - BonA);

end



%% determine amplitude adjusted pressure and flow at each segment
function [A_adjusted] = adjust_amplitude(A, BonA)
% The first segment sends an outgoing wave with amplitude A = 1;
% the reflection coefficent R = B/A means B/A is reflected at each segment
% so the next segment has an input wave amplitude A* = 1-BonA
A_adjusted = A.*BonA;

end

%% General transfer function calculation %%%%%%%%% not yet used %%%%%%%%%%
function [TF] = transfer_fn(specified_output, specified_input)
% care should be taken to choose appropriate real/abs values
TF = specified_output/specified_input;

end

%% BANDWIDTH MEASUREMENTS (Optional)
function [output] = formant_bandwidth(PARAMS, Transfer)

freq = PARAMS.freq;
resolution = PARAMS.resolution;

% need to specify upper and lower bounds of the frequency
% for each pair of max and min
max_first = [150 1050 2250 3250];
max_last = [1000 2200 3200 4000];
% create vectors of the right size to be filled below
z_max = [nan nan nan nan];
f_max = [nan nan nan nan];
max_bandwidth = [nan nan nan nan];

for n = 1:4
    % convert to index of the frequency vector
    [~, c2] =find(freq>=(max_first(n)));
    low2(n) = min(c2); % max_first index above the extreme
    [~, c2] =find(freq<=(max_last(n)));
    high2(n) = max(c2); % max_first index above the extreme
    
    [value2, index2] = max(abs(Transfer(low2(n):high2(n))));
    % zmax is the maximum impedance value
    z_max(n) = value2;
    % fmax is the frequency of the minimum impedance
    f_max(n) = max_first(n) - resolution + (index2.*resolution);
    [~, max_index] = find(abs(Transfer(low2(n):high2(n)))...
        >=(max(abs(Transfer(low2(n):high2(n))))./sqrt(2)));
    f_max_range = max_first(n) - resolution + (max_index.*resolution);
    if max(f_max_range) >= 1
        max_bandwidth(n) = max(f_max_range) - min(f_max_range);
    end
    
    % ignore max_bandwidth greater than 250
    if max_bandwidth(n) > 250;
        max_bandwidth(n) = nan;
    end
    
end

% OUTPUT
output.T = Transfer;
output.max.z(:,1) = z_max;
output.max.f(:,1) = f_max;
output.max.b(:,1) = max_bandwidth;

end



%% BANDWIDTH MEASUREMENTS (Optional)
function [output] = impedance_bandwidth(PARAMS, Transfer)

freq = PARAMS.freq;
resolution = PARAMS.resolution;

% need to specify upper and lower bounds of the frequency
% for each pair of max and min
max_first = [150 700 1800 2600 3600];
max_last = [400 1500 2500 3500 4200];
% create vectors of the right size to be filled below
z_max = [nan nan nan nan nan];
f_max = [nan nan nan nan nan];
max_bandwidth = [nan nan nan nan nan];

% min_first = [10 200 1000 2000 3100];
% min_last = [100 800 2000 3000 4000];

min_first = [10 400 1100 2000 3000 ];
min_last = [100 800 1800 2600 3600 ];

% create vectors of the right size to be filled below
z_min = [nan nan nan nan nan];
f_min = [nan nan nan nan nan];
min_bandwidth = [nan nan nan nan nan];

for n = 2:5 %%%%%%%%%%%%%% change this to look at other resonances
    % convert to index of the frequency vector
    [~, c2] =find(freq>=(max_first(n)));
    low2(n) = min(c2); % max_first index above the extreme
    [~, c2] =find(freq<=(max_last(n)));
    high2(n) = max(c2); % max_first index above the extreme
    
    [value2, index2] = max(abs(Transfer(low2(n):high2(n))));
    % zmax is the maximum impedance value
    z_max(n) = value2;
    % fmax is the frequency of the minimum impedance
    f_max(n) = max_first(n) - resolution + (index2.*resolution);
    [~, max_index] = find(abs(Transfer(low2(n):high2(n)))...
        >=(max(abs(Transfer(low2(n):high2(n))))./sqrt(2)));
    f_max_range = max_first(n) - resolution + (max_index.*resolution);
    if max(f_max_range) >= 1
        max_bandwidth(n) = max(f_max_range) - min(f_max_range);
    end
    
    % % %     % limit max_bandwidth to 250
    % % %     if max_bandwidth(n) > 250;
    % % %         max_bandwidth(n) = nan;
    % % %     end
    
    
    % convert to index of the frequency vector
    [~, c1] =find(freq>=(min_first(n)));
    low1(n) = min(c1); % min_first index above the extreme
    [~, c1] =find(freq<=(min_last(n)));
    high1(n) = max(c1); % min_first index above the extreme
    
    
    [value1, index1] = min(abs(Transfer(low1(n):high1(n))));
    % zmin is the minimum impedance value
    z_min(n) = value1;
    % fmin is the frequency of the minimum impedance
    f_min(n) = min_first(n) - resolution + (index1.*resolution);
    [~, min_index] = find(abs(Transfer(low1(n):high1(n)))...
        <=(min(abs(Transfer(low1(n):high1(n)))).*sqrt(2)));
    f_min_range = min_first(n) - resolution + (min_index.*resolution);
    if min(f_min_range) >= 1
        min_bandwidth(n) = max(f_min_range) - min(f_min_range);
    end
    
    % % %     % limit the min_bandwidth to 250
    % % %     if min_bandwidth(n) > 250;
    % % %         min_bandwidth(n) = nan;
    % % %     end
    
end

% OUTPUT
% Choose which variables to save to the output
output.T = Transfer;
output.max.z(:,1) = z_max;
output.max.f(:,1) = f_max;
output.max.b(:,1) = max_bandwidth;

output.min.z(:,1) = z_min;
output.min.f(:,1) = f_min;
output.min.b(:,1) = min_bandwidth;

end



