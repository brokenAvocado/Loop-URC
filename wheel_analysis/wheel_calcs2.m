clear;
clc;

x_R = 1.2;
z_R = 1.2;
y_R = 0.8;

% If using Markforged MK2 3 Printer max volume (12.6 x 5.2 x 6 in)

% Initialize values for global variables
m = 25; % mass/0.7 for weight margin analysis (AKA design rover to be 70% of this mass value) # rover mass {kg}
b = 4/39.37; % width of wheels {m}
D = 8/39.37; % wheel diameter {m}
N = 6; % qty of wheels
w_s = 0.275; % wheel spacing from edge of grouser to edge of grouser {m}
g = 9.81; % gravity on earth {m/s^2}

% Grousers
y_g = 0.006; % height of grouser cross sect face {m} (i.e. looking at grouser head-on from the part that touches the ground, this is the height of its face)
N_Grouser = 60; % total number of grousers

% Dynamics
v = 0.2; % rover speed {m/s}

SF = 1.4; % see NASA-STD-5001B, section: Metllatic Structures (page 18)

%% Material Properties

% % Makerforged® Onyx™ (Carbon Fiber Composite) 3D Printing Filament (https://www-objects.markforged.com/craft/materials/CompositesV5.2.pdf)
% Su = 37e6; % Ultiamte tensile strength {Pa}
% Sy = 40e6; % Yield Strength {Pa}
% S_fatigue = 27.75e6; % Fatigue strength~ Estimate failure at 75% {Pa} based of: https://www.sciencedirect.com/science/article/pii/S1359836821000974#sec3 and https://www.sciencedirect.com/science/article/pii/S1359836821000974#sec3
% E = 2.4e9; % Modulus of Elasticity {Pa}
% G = 3e9; % Shear modulus
% Poisson =  0.5*(E/G)-1; % Poisson's ratio
% Density = 1200; % kg/m^3

% Ultimaker 95A TPU (https://3dneworld.com/wp-content/uploads/2017/08/TPU95A.pdf)
Su = 39e6; % Ultiamte tensile strength {Pa}
Sy = 8.6e6; % Yield Strength {Pa}
S_fatigue = 4.3e6; % Fatigue strength~ Estimate failure at 75% {Pa} based of: https://www.sciencedirect.com/science/article/pii/S1359836821000974#sec3 and https://www.sciencedirect.com/science/article/pii/S1359836821000974#sec3
E = 26e6; % Modulus of Elasticity {Pa}
G = 0.102e9; % Shear modulus (https://www.matweb.com/search/datasheettext.aspx?matguid=9f5318a1f93b403bbd5748abec70fac1#:~:text=Shear%20Modulus%2C%20Average%20value%3A%200.102%20GPa%20Grade%20Count%3A3)
Poisson =  0.5*(E/G)-1; % Poisson's ratio
Density = 1200; % kg/m^3

%% Soil Characteristics

k_c = 74600; % modulus of cohesion of soil deformation {N/m^(n+1)}
k_phi = 2080000; % modulus of friction of soil deformation {N/m^(n+2)}
k = (k_c / b) + k_phi; % modulus of soil deformation

n = 1.1;
phi = 33.7; % angle of internal resistace of soil {deg} 
c = 3300; % Cohesion {N/m^2}

N_q = exp((((3*pi)/2)-deg2rad(phi))*tand(phi))/(2*(cosd(rad2deg(pi/4) + (phi)/2))^2);
N_c = (N_q-1)/tand(phi);
N_gamma = (2*(N_q+1)*tand(phi))/(1+(0.4*sind(4*phi)));

gamma = 2315.625; % weight density of soil {N/m^3}

K_bar_c = (N_c - tand(phi))*((cosd(phi))^2); % Cohesive modulus of soil deformation 
K_bar_phi = 72.77; % Frictional modulus of soil deformation
K_bar_gamma = (((2*N_gamma) / (tand(phi))) + 1) * ((cosd(phi))^2);
K_bar = 0.018; %shear deformation modulus {m}

theta = 25; % maximum slope angle the rover is expected to traverse upwards {deg}

%% Parameters

W = m*g; % gross weight {N}
W_w = W / N; % weight on wheel {N}

z = ((3*W_w) / ((3-n)*b*k*(D^.5)))^(2/(2*n+1)); % distance below ground surface (i.e. from surface to bottom of wheel) {m}
l = (D/2) * acos(1 - (2*z)/D); % length of contact patch {m}
A = b * l; %contact area of wheel on soil {m^2} 

h = (z_R - (w_s * ((N / 2) - 1)) - (D*(N/2))) / N; % derive grouser height from wheel diameter and maximum size possible based on size constraint and wheel count
x_g = b; % width of grouser cross sect face {m} (i.e. looking at grouser head-on from the part that touches the ground, this is the width of its face) (equal to wheel width)
N_g = (N_Grouser*l) / ((2*pi)*(D/2)); % number of grousers in contact with the ground (based on arc length of the wheel in contact with the ground and the count of grousers)

% round down if decimal is not at least 0.8 or more
% prety arbitrary but conservative (want 80% contact to say the grouser is in contact with the soil)
if (N_g - floor(N_g)) >= 0.8
    N_g = round(N_g); %round to nearest whole number
else
    N_g = floor(N_g); %round to floor
end

alpha = acosd(1 - ((2*z)/D)); % angle of attack of wheel in soil {deg}

l_o = z*(tand(rad2deg(pi / 4)) - phi/2)^2;

s = 0.5; % s = 1 - v / (omega * (D/2)) % wheel slip ratio (assume a conservative value of 0.5)
omega = v / ((1 - s)*(D/2)); % omega = (v / (D/2)) * 1.05 % solve for omega from wheel slip ratio

length = (D*(N/2) + 2*h*(N/2) + w_s * ((N / 2) - 1));

%% Constraint Check

if (length <= z_R)
    fprintf("PASS: Wheel Diameter, Wheel Count, and Grouser Height are within allowances. %fm is <= z_R = %fm\n", length, z_R);
else
    fprintf("FAIL: Wheel Diameter, Wheel Count, and Grouser Height FAIL the contraint check. One of these is too large. %fm is > z_R = %fm\n", length, z_R);
end

% weight constraint
W_max_per_wheel = A*(c*N_c + gamma*z*N_q + (1/2)*gamma*b*N_gamma); % max weight on soil per wheel {n}
W_max = N*W_max_per_wheel;

if (W < W_max)
    fprintf("PASS: Rover weight and wheel thickness are within allowances. Weight: %f < W_max: %f\n", W, W_max);
else
    fprintf("FAIL: Rover weight and wheel thickness are NOT within allowances. Weight may be too high, or thickness is too small. Weight: %f >= W_max: %f\n", W, W_max);
end

% grouser cross sectional y_g height or "thickness" constraint 
distance_between_grousers = (2*pi*(D/2))/N_Grouser; % spacing between grousers {m}

if (y_g < distance_between_grousers)
    fprintf("PASS: Grouser y_g cross sectional face height is within allowances. y_g = %fm is < %fm\n", y_g, distance_between_grousers);
else
    fprintf("FAIL: Grouser y_g cross sectional face height is NOT within allowances. y_g = %fm is >= %fm\n", y_g, distance_between_grousers);
end

% grousers in contact constraint
if (N_g >= 1)
    fprintf("PASS: Number of Grousers in contact with soil: %d is >= 1\n", N_g);
else
    fprintf("FAIL: Number of Grousers in contact with soil: %d is < 1\n", N_g);
end

% ignore bulldozing resistance constraint
% based on apostolopoulos_dimitrios_2001_dissertation
ratio = z/D;
if (ratio < 0.06)
    fprintf("PASS: Can ignore Bulldozing Resistance: %f is < 0.06\n", ratio);
else
    fprintf("FAIL: CANNOT ignore Bulldozing Resistance: %f is >= 0.06\n", ratio);
end

fprintf("_______________________________________ \n")

R_c = compressionResistance(b, k, n, z);
fprintf("Compression Resistance = " + R_c + " N\n")

R_r = rollingResistance(b, k, W, D);
fprintf("Rolling Resistance = " + R_r + " N\n")

R_g = gravitationalResistance(W, theta);
fprintf("Gravitational Resistance = " + R_g + " N\n")

R_b = bullzodeResistance(alpha, phi, b, z, c, K_bar_c, gamma, K_bar_gamma, l_o, D);
fprintf("Bullzode Resistance = " + R_b + " N\n")

H_smooth = smoothTractiveForce(A, c, W_w, phi, K_bar, s, l, N);
fprintf("Smooth Tractive Force = " + H_smooth + " N\n")

H_grouser = grouserTractiveForce(b, l, c, h, N_g, W_w, phi, K_bar, s, N);
fprintf("Grouser Tractive Force = " + H_grouser + " N\n")

DP_smooth = drawbarPullSmooth(H_smooth, R_c, R_b, R_g, R_r);
fprintf("Drawbar Pull Smooth = " + DP_smooth + " N\n")

DP_grouser = drawbarPullGrouser(H_grouser, R_c, R_b, R_g, R_r);
fprintf("Drawbar Pull Grouser = " + DP_grouser + " N\n")

torqueReq = torqueReqGrouser(H_grouser, N, D);
fprintf("Torque Required = " + torqueReq + " N-m\n")

powerReq = powerReqPerWheel(torqueReq, omega);
fprintf("Power Required = " + powerReq + " W\n")

fprintf("_______________________________________ \n")

%% Stress on Grousers

Force1 = (H_grouser/N)/N_g; % Thrust load {N}
distance = h/2; % h/2 since this is a distributed load
M1 = Force1*distance; % force times distance (height of grouser times drawbar pull of grouser / number of grousers in contact)
Ix = (1/12)*(x_g)*(y_g^3); % moment of intertia of grouser in x {m^4}
c_g1 = y_g / 2; %height from neutral axis to load {m}
Stress_b1 = M1*c_g1 / Ix;

Force2 = 40.9/N_g; % Side load  {N} (40.9 is placeholder from steering Eqn thats in progress)
M2 = Force2*distance; % force times distance (height of grouser times drawbar pull of grouser / number of grousers in contact)
Iy = (1/12)*(y_g)*(x_g^3); % moment of intertia of grouser in x {m^4}
c_g2 = x_g / 2; % height from neutral axis to load {m}
Stress_b2 = M2*c_g2 / Iy;

Stress_b = sqrt(Stress_b1^2 + Stress_b2^2);

MOS_yield = (Sy / (Stress_b*SF)) - 1;
fprintf("Margin of Safety wrt Yield Due to Bending: %f\n", MOS_yield);

MOS_fatigue = (S_fatigue / (Stress_b*SF)) - 1;
fprintf("Margin of Safety wrt fatigue Due to Bending: %f\n", MOS_fatigue);

V_force1 = (H_grouser/N)/N_g; % Thrust load {N}
Stress_shear1 = V_force1 / (y_g*x_g); % Force per grouser / cross sect area

V_force2 = 9.92/N_g; % Side force {N} (given by steering in the section below)
Stress_shear2 = V_force2 / (y_g*x_g); % Force per grouser / cross sect area

% Calculation of combined loading
Stress_shear = sqrt(Stress_shear1^2 + Stress_shear2^2);

% Calculation of margin of safety
MOS_yield = (Sy / (Stress_shear*SF)) - 1;
fprintf('Margin of Safety wrt Yield Due to Shear: %f\n', MOS_yield);

MOS_fatigue = (S_fatigue / (Stress_shear*SF)) - 1;
fprintf('Margin of Safety wrt fatigue Due to Shear: %f\n', MOS_fatigue);

K_euler = 2; % fixed-free
Pcr = (pi * E * Ix) / (K_euler*h)^2;

% buckling load is greater force between the applied tractive force, H_grouser, divided by number of grousers in contact and the weight of the rover divided by the number of grousers in contact
% This assumes that the worse case condition is either entirely tractive force causes buckling (in the case of a grouser that is paralell wrt the soil buckling against a soil wall)
% or the weight of the rover on the grouser as it is direction perpindular wrt the soil

P_load = (36.7/2)/N_g; % Suspension load down worst incline {N} (given by suspension in the section below)

if ((H_grouser/N_g) > P_load)
    P_load = H_grouser/N_g;
end

MOS_buckling = (Pcr / (P_load*SF)) - 1;
fprintf("Margin of Safety wrt Buckling: " + MOS_buckling + "\n");

Deflection_grouser1 = Force1*(h^3) / (3*E*Ix); %PL^3 / (3EI) {m}
Deflection_grouser2 = Force2*(h^3) / (3*E*Iy); %PL^3 / (3EI) {m}

%find the greatest deflection
if (Deflection_grouser1 < Deflection_grouser2)
Deflection_grouser = Deflection_grouser2;
else
Deflection_grouser = Deflection_grouser1;
end

fprintf("Deflection of grouser in m: " + Deflection_grouser + "\n");

% Data_amount = 50;
% X = linspace(0.01,0.5,Data_amount);
% 
% for k = 1:Data_amount
%     b=X(k);
%     R_c = compressionResistance(b, k, n, z);
%     R_r = rollingResistance(b, k, W, D);
%     R_g = gravitationalResistance(W, theta);
%     R_b = bullzodeResistance(alpha, phi, b, z, c, K_bar_c, gamma, K_bar_gamma, l_o, D);
%     H_grouser = grouserTractiveForce(b, l, c, h, N_g, W_w, phi, K_bar, s, N);
%     
%     y(k) = drawbarPullGrouser(H_grouser, R_c, R_b, R_g, R_r);
% end
% plot(X, y);

%% Functions for Wheels

function R_c = compressionResistance(b,k,n,z)
  R_c_one = (b*k/(n+1))*(z^(n+1)) ;

  %compression resistance for entire rover (2 sides to a rover)
  R_c = R_c_one * 2 ;

end

function R_r = rollingResistance(b,k,W,D)
   R_r = (0.580 / ((b*k)^.5)) * ((W^(3/2)) / (D^(3/4)));
end

function R_g = gravitationalResistance(W,theta)
   R_g = W*sind(theta);
end

function R_b = bullzodeResistance(alpha, phi, b, z, c, K_bar_c, gamma, K_bar_gamma, l_o, D)
  % Compression resistance per leading wheel {N}
  R_b_leading = ((b * sind(alpha + phi)) / (2 * sind(alpha) * cosd(phi))) * (2 * z * c * K_bar_c + gamma * (z^2) * K_bar_gamma) + ((l_o^3) * (gamma / 3)) * ((pi/2) - deg2rad(phi)) + c * (l_o^2) * (1 + tand((rad2deg(pi/4) + (phi/2))));

  % Total bulldozing resistance (assume 2 leading wheels)
  R_b = R_b_leading * 2; % {N}
  
  if (z/D < 0.06)
    R_b = 0;
  end
end

function H_smooth = smoothTractiveForce(A,c,W_w,phi,K_bar,s,l,N)
H_smooth_per_wheel = (A*c + W_w*tand(phi)) * (1 - (K_bar/(s*l))*(1 - exp((-1*s*l)/K_bar)));

H_smooth = N * H_smooth_per_wheel;
end

function H_grouser = grouserTractiveForce(b,l,c,h,N_g,W_w,phi,K_bar,s,N)
  H_grouser_per_wheel = (b*l*c*(1+2*h/b)*N_g + W_w*tand(phi)*(1 + 0.64*(h/b)*atan(b/h)))*(1 - (K_bar / (s*l))*(1 - exp((-s*l)/K_bar)));

  H_grouser = N * H_grouser_per_wheel;
end

function DP = drawbarPullSmooth(H_smooth, R_c, R_b, R_g, R_r)
  DP = H_smooth - (R_c + R_b + R_g + R_r);
end

function DP = drawbarPullGrouser(H_grouser, R_c, R_b, R_g, R_r)
  DP = H_grouser - (R_c + R_b + R_g + R_r);
end

function torqueReq = torqueReqGrouser(H_grouser, N, D)
  torqueReq = (H_grouser/N)*(D/2); %Force divided by qty of wheels times distance
end

function powerReq = powerReqPerWheel(torque,omega)
  powerReq = torque * omega; % power times rpm
end