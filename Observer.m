
clc
clear all
close all

load GNSSaidedINS_data.mat
rng('default');
%% System Inputs initialization
% v=3;
% omega=[1;2;3];
%% Parameters Initialization
numt = 299; % time steps of the trajectory
tspan=0.01; % time period of integration
d=7; % Dimensionality of the state vector
muBar_state=zeros(7,1,numt); % Mean state Vector
real_obv=zeros(6,numt); % Measurements from the real system
state_trajec=zeros(7,numt); % State Measurements from the real system
observation=zeros(6,1,2*d+1,numt); % Observation vectors for each sigma point
sv_sp=zeros(7,1,2*d+1); % State vectors for each sigma point
sp=zeros(7,1,2*d+1); % Sigma Points
sp_tr=zeros(7,1,2*d+1); % Sigma Points after the first transformation
muBar_obv=zeros(6,1,numt); % Mean observation vector
InoCov=zeros(6,6); % Innovation Covariance Martrix
SO_Cov=zeros(7,6); % State-Observation Covariance matrix
mu_t=zeros(7,numt); % Mean State vector matrix
a= zeros(7,1,1,1); % Buffer variable for Mean of State vector Sigma Points
b= zeros(7,7);     % Buffer variable for Covariance of the transformed Sigma Points
q= zeros(6,1,1,1); % Buffer variable for Mean Observation
L=zeros(6,6,1,1);  % Buffer variable for Innovation covarianc
x1=zeros(7,6);     % Buffer variable for State-Observation Covariance
R=0.01*eye(7); % Noise covariance of process
Q=diag([.001,.001,.001,.001,.001,.001]);% Noise covariance of observation (measurement)
 % Actual initial state vector

velocity = in_data.SPEEDOMETER.speed(1,1:numt); %get 300 values of speed
time=in_data.SPEEDOMETER.t(1:300,:);% define time for 300 instances
time = time';
%% 
%get acceleration data from velocity
acc = zeros(1,numt);
for i=1:numt-1
    acc(i)=(velocity(i+1)-velocity(i))/(time(i+1)-time(i));
end

%%
%get transversal and forward accelerometer reading
trans_x = in_data.IMU.acc(1,1:numt);
frwd_y = in_data.IMU.acc(2,1:numt);
angvel_z = in_data.IMU.gyro(3,1:numt);

g = 9.81;
psi(1,1) = 11.02;
lamda(1,1) = 77.00;
alt(1,1) = 0;
vel(1,1) = 0;
pit(1,1) = asin(frwd_y(1,1) - acc(1,1)/g);
rol(1,1) = -asin(trans_x(1,1)) + velocity(1,1) * angvel_z(1,1)/g;
azi(1,1) = atan2(trans_x(1,1),frwd_y(1,1));
omega_e = 7.2921159e-5; %earth rotation rate
R_M = 6378137;
R_N = 6356752.3;

mu_treal=[psi(1,1);lamda(1,1);alt(1,1);vel(1,1);pit(1,1);rol(1,1);azi(1,1)];

%%  trajectory Generator

for i=1:299
    velo_noise(:,i) = velocity(:,i) + (rand(1,1));
    acc_noise(:,i) = acc(:,i) + (rand(1,1));
    trans_x_noise(:,i) = trans_x(:,i) + (rand(1,1));
    frwd_y_noise(:,i) = frwd_y(:,i) + (rand(1,1));
    angvel_z_noise(:,i) = angvel_z(:,i) + (rand(1,1));
end

%%
%form state transition matrix

for k = 2:numt
    velo_f(:,k) = velocity(:,k-1) - velo_noise(:,k-1);
    format long;
    psi(:,k)= psi(:,k-1) + (velo_f(:,k-1).*cos(azi(:,k-1)).*cos(pit(:,k-1))./ R_M + alt(:,k-1)).*tspan;
    lamda(:,k) = lamda(:,k-1) + (velo_f(:,k-1).*sin(azi(:,k-1)).*cos(pit(:,k-1))./(R_N + alt(:,k-1).*cos(psi(:,k-1)))).*tspan;
    alt(:,k) = alt(:,k-1) + (velo_f(:,k-1).*sin(pit(:,k-1))).*tspan;
    vel(:,k) = velocity(:,k-1) - velo_noise(:,k-1);
    pit(:,k) = asind((trans_x(:,k-1) - trans_x_noise(:,k-1))-(acc(:,k-1)-acc_noise(:,k-1)))./g;
    rol(:,k) = -asind((frwd_y(:,k-1) - frwd_y_noise(:,k-1)) + (velocity(:,k-1) - velo_noise(:,k-1)).*(angvel_z(:,k-1) - angvel_z_noise(:,k-1)))./g.*cos(pit(:,k));
    gamma_z(:,k) = (angvel_z(:,k-1) - angvel_z_noise(:,k-1));
    U_E(:,k) = sin(azi(:,k-1)).*cos(pit(:,k-1)).*cos(gamma_z(:,k-1)) - (cos(azi(:,k-1)).*cos(rol(:,k-1))+sin(azi(:,k-1)).*sin(pit(:,k-1)).*sin(rol(:,k-1))).*sin(gamma_z(:,k-1)).*tspan;
    U_N(:,k) = cos(azi(:,k-1)).*cos(pit(:,k-1)).*cos(gamma_z(:,k-1)) - (-sin(azi(:,k-1)).*cos(rol(:,k-1))+cos(azi(:,k-1)).*sin(pit(:,k-1)).*sin(rol(:,k-1))).*sin(gamma_z(:,k-1)).*tspan;
    azi(:,k) = atan(U_E(:,k)/U_N(:,k)) + omega_e.*sin(psi(:,k-1)).*tspan+ (velo_f(:,k-1).*sin(azi(:,k-1)).*cos(pit(:,k-1)).*tan(psi(:,k-1))./R_N + alt(:,k-1)).*tspan;
end

x_k = [psi;lamda;alt;vel;pit;rol;azi];

%%rotation matrix
rot_mat = [cos(azi(:,k-1)).*cos(rol(:,k-1))+sin(azi(:,k-1)).*sin(pit(:,k-1)).*sin(rol(:,k-1))   sin(azi(:,k-1)).*cos(pit(:,k-1))   cos(azi(:,k-1)).*sin(rol(:,k-1))-sin(azi(:,k-1)).*sin(pit(:,k-1)).*cos(rol(:,k-1));
    -sin(azi(:,k-1)).*cos(rol(:,k-1))+cos(azi(:,k-1)).*sin(pit(:,k-1)).*sin(rol(:,k-1))       cos(azi(:,k-1)).*cos(pit(:,k-1))     -sin(azi(:,k-1)).*sin(rol(:,k-1))+cos(azi(:,k-1)).*sin(pit(:,k-1)).*cos(rol(:,k-1));
    -cos(pit(:,k-1)).*sin(rol(:,k-1))                                                           sin(pit(:,k-1))                     cos(pit(:,k-1)).*cos(rol(:,k-1))];


%% weights and Covariance Initialization
[wt_pt,wt_cov]=sigmaPntCovWt(d);
covariance(1:7,1:7)= eye(7);
mu_t(1:7,1) = x_k(1:7,1)+ 0.001*randn(7,1);% try removing randn and simulate
%% Filter
for t=2:numt %% iterations for time steps equal to 0.025 seconds
    %% Sigma Points for each state vector
    covariance(1:7,1:7)= 0.5*eye(7);
    sp(1:7,1,1:2*d+1)=sigmapnt(mu_t(1:7,t-1),covariance(1:7,1:7),d);
    
    %% Propogate sigma points into the process model
    
    for k=2:2*d+1 %%iteration for each sigma point
     velo_f_sigma = zeros(size(sp(1,:,k)));
     %psi_sigma = zeros(size(sp(1,:,k)));
    % pit_sigma = zeros(size(sp(1,:,k)));
    velo_f_sigma(1,:,k) = sp(4,:,k);
    format long;
    psi_sigma(1,:,k-1)= sp(1,:,k-1) + (velo_f_sigma(1,:,k-1).*cos(sp(7,:,k-1)).*cos(sp(5,:,k-1))./ R_M + sp(3,:,k-1)).*tspan;
    lamda_sigma(1,:,k-1) = sp(2,:,k-1) + (velo_f_sigma(1,:,k-1).*sin(sp(7,:,k-1)).*cos(sp(5,:,k-1)))./(R_N + sp(3,:,k-1)).*cos(sp(1,:,k-1)).*tspan;
    alt_sigma(1,:,k-1) = sp(3,:,k-1) + (velo_f_sigma(1,:,k-1).*sin(sp(5,:,k-1))).*tspan;
    vel_sigma(1,:,k-1) = sp(4,:,k-1);
    pit_sigma(1,:,k-1) = asind((trans_x(:,k-1) - trans_x_noise(:,k-1))-(acc(:,k-1)-acc_noise(:,k-1)))./g;
    rol_sigma(1,:,k-1) = -asind((frwd_y(:,k-1) - frwd_y_noise(:,k-1)) + (velocity(:,k-1) - velo_noise(:,k-1)).*(angvel_z(:,k-1) - angvel_z_noise(:,k-1)))./g.*cos(sp(5,:,k));
    gamma_z(:,k) = (angvel_z(:,k-1) - angvel_z_noise(:,k-1));
    U_E(:,k) = sin(sp(7,:,k-1)).*cos(sp(5,:,k-1)).*cos(gamma_z(:,k-1)) - (cos(sp(7,:,k-1)).*cos(sp(6,:,k-1))+sin(sp(7,:,k-1)).*sin(sp(5,:,k-1)).*sin(sp(6,:,k-1))).*sin(gamma_z(:,k-1)).*tspan;
    U_N(:,k) = cos(sp(7,:,k-1)).*cos(sp(5,:,k-1)).*cos(gamma_z(:,k-1)) - (-sin(sp(7,:,k-1)).*cos(sp(6,:,k-1))+cos(sp(7,:,k-1)).*sin(sp(5,:,k-1)).*sin(sp(6,:,k-1))).*sin(gamma_z(:,k-1)).*tspan;
    azi_sigma(1,:,k-1) = atan(U_E(:,k)/U_N(:,k)) + omega_e.*sin(sp(1,:,k-1)).*tspan+ (velo_f(:,k-1).*sin(sp(7,:,k-1)).*cos(sp(5,:,k-1)).*tan(sp(1,:,k-1))./R_N + sp(3,:,k-1)).*tspan;      
    sv_sp=[psi_sigma;lamda_sigma;alt_sigma;vel_sigma;pit_sigma;rol_sigma;azi_sigma];
    end
%sv_sp(1,:,:) = sp(1,:,1:14);

    %% Mean of State vector Sigma Points
    a = zeros(7,1,1,1);
    for msgmp=1:2*d
%         a = a+wt_pt(msgmp)*sv_sp(1:7,1,msgmp);
    a = a+wt_pt(msgmp)*sv_sp(1:7,1,msgmp);
    end
    muBar_state(:,1,t) = a;
%    x_k_mean = zeros(7,1,1,1);
% [x_k_lat,x_k_lon] = meanm(x_k(1,:),x_k(2,:));
% x_k_alt = mean(x_k(3,:))./length(t+1);
% x_k_vel = mean(x_k(4,:))./length(t+1);
% x_k_pit = mean(x_k(5,:))./length(t+1);
% x_k_rol = mean(x_k(6,:))./length(t+1);
% x_k_azi = mean(x_k(7,:))./length(t+1);
% 
% x_k_mean = [x_k_lat;x_k_lon;x_k_alt;x_k_vel;x_k_pit;x_k_rol;x_k_azi];
% muBar_state(:,1,t) = x_k_mean;

    %% Covariance of the transformed Sigma Points
    b= zeros(7,7);
    for msgmcov=1:2*d
        covariance(1:7,1:7)= b+wt_cov(msgmcov)*(sv_sp(1:7,1,msgmcov)-muBar_state(1:7,1,t))*(sv_sp(1:7,1,msgmcov)-muBar_state(1:7,1,t))';
        b=covariance(1:7,1:7);
    end
    covariance(:,:)=covariance(:,:);
    covariance = 0.5*eye(7,7);
    %% Sigma points around the transformed mean
    
    sp_tr(1:7,1,1:2*d+1)=sigmapnt(muBar_state(1:7,1,t),covariance(1:7,1:7),d);
    
    %% Observation matrix determination
earth = referenceSphere('Earth');
lat0=11.0243584;
long0=77.0061533;
h0=0;

xNorth=in_data.GNSS.pos_ned(1,:);
xEast=in_data.GNSS.pos_ned(2,:);
xDown=in_data.GNSS.pos_ned(3,:);
[lat,lon,h]=ned2geodetic(xNorth,xEast,xDown,lat0,long0,h0,earth); 

   for k=1:299
    z_lat(:,k) = lat(:,k) + rand(1,1)*0.001;
    z_lon(:,k) = lon(:,k) + rand(1,1)*0.001;
    z_alt(:,k) = h(:,k) + rand(1,1)*0.001;
    zvel_east(:,k) = (velo_f(:,k)).*sin(azi(:,k)).*cos(pit(:,k)) + rand(1,1);
    zvel_north(:,k) = (velo_f(:,k)).*cos(azi(:,k)).*cos(pit(:,k)) + rand(1,1);
    zvel_up(:,k) = (velo_f(:,k)).*sin(pit(:,k)) + rand(1,1);
   end
z_k = [z_lat;z_lon;z_alt;zvel_east;zvel_north;zvel_up];
for i = 1:299
    ned_vel = [zvel_east;zvel_north;zvel_up];
    vel_body(:,i) = inv(rot_mat)*ned_vel(:,i);
end
    %% Mean Observation
    z_cov = eye(6,6);
    e = 6;
    observation(1:6,1,1:2*e+1)=sigmapnt(z_k(1:6,t-1),z_cov(1:6,1:6),e);
 %%
%     for k=1:13
%     z_lat_sigma(:,1,k,t) = sp_tr(1,:,k) ;
%     z_lon_sigma(:,1,k,t) = sp_tr(2,:,k) ;
%     z_alt_sigma(:,1,k,t) = sp_tr(3,:,k) ;
%     zvel_east_sigma(:,1,k,t) = (velo_f(:,k)).*sin(sp_tr(7,:,k)).*cos(sp_tr(5,:,k));
%     zvel_north_sigma(:,1,k,t) = (velo_f(:,k)).*cos(sp_tr(7,:,k)).*cos(sp_tr(5,:,k)) ;
%     zvel_up_sigma(:,1,k,t) = (velo_f(:,k)).*sin(sp_tr(5,:,k)) ;
%     end
%     observation = [z_lat_sigma;z_lon_sigma;z_alt_sigma;zvel_east_sigma;zvel_north_sigma;zvel_up_sigma];
    q=zeros(6,1,1,1);
    for msgm_obv=1:2*e+1
        q=q+wt_pt(msgm_obv)*observation(1:6,1,msgm_obv);
    end
    muBar_obv(:,1,t)=q;
    %% Innovation covariance
    L=zeros(6,6,1,1);
    for f=1:2*e+1
        L=L+ wt_cov(f)*(observation(:,1,f)-muBar_obv(:,1,t))*(observation(:,1,f)-muBar_obv(:,1,t))';
    end
    InoCov(:,:)=L;
    InoCov(:,:)=InoCov(:,:) ;
   InoCov = 1.225*eye(6,6);
    %% State-Observation Covariance
    x1=zeros(7,6);
    for g=1:2*d+1
        x1= x1+ wt_cov(g)*(sp_tr(:,1,g)-muBar_state(:,1,t))*(observation(:,1,g)-muBar_obv(:,1,t))';
    end
    SO_Cov(:,:)=x1;
  SO_Cov(:,:)= 1.225*eye(7,6);
    %% Kalman Gain
    K(1:7,1:6,t)=SO_Cov(:,:)*(InoCov(:,:))^-1;
   % real_obv=zeros(6,numt); % Measurements from the real system
    %% Mean State
       mu_t(:,t)=muBar_state(:,1,t)+K(:,:,t)*(real_obv(:,t)-muBar_obv(:,1,t));
   
    %% Covariance SP of the states
    covariance(:,:)=covariance(1:7,1:7)-K(:,:,t)*InoCov(:,:)*K(:,:,t)';
end





