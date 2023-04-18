%% Drive Wheel Analysis, April 16, 2023

% Developer Notes: K, soil slip, is unknown, assuming lunar soil (1.78
% cub cm)

%% Constants

% Bounding Box Dimensions (all in meters)
x_R = 1.2;
z_R = 1.2;
y_R = 0.8;

SF = 1.4; % see NASA-STD-5001B, section: Metllatic Structures (page 18)

%% Soil Characteristics
% We chose a sandy loam #2 on this table, https://tinyurl.com/ycydbntc

% Shibly's soil parameters
n = 0.7; % some non dimensional factor
k_c = 5270; % modulus of cohesion of soil {N/m^(n+1)}
k_phi = 1515040; % modulus of friction of soil {N/m^(n+2)}
k = (k_c / b) + k_phi; % modulus of soil deformation
phi = 33; % angle of internal resistace of soil {deg}
c = 170; % cohesion {N/m^2}

% Terzaghi's soil theory
N_q = exp(2*pi*(0.75-phi/360)*tand(phi))/(2*cosd(45+phi/2)^2); % 32.23
N_c = (N_q-1)/tand(phi);
N_gamma = 2*(N_q+1)*tand(phi)/(1+0.4*sind(4*phi));

gamma = 15200.31; % weight density of soil {N/m^3}

% Bulldozing parameters (unitless and based on Terzaghi)
K_bar_c = (N_c-tand(phi))*cosd(phi)^2; % cohesive modulus of soil deformation 
K_bar_gamma = (((2*N_gamma) / (tand(phi))) + 1) * (cosd(phi))^2; % density modulus of soil deformation
K_bar = 0.018; % coefficient of soil slip {m}

theta = 25; % maximum uphill angle {deg}

%% Parameters
W = m*g; % gross weight {N}
W_w = W / N; % weight on wheel {N}

z = ((3*W_w) / ((3-n)*b*k*(D^.5)))^(2/(2*n+1)); % distance below ground surface (i.e. from surface to bottom of wheel) {m}
l = (D/2) * acos(1 - (2*z)/D); % length of contact patch {m}
A = b * l; %contact area of wheel on soil {m^2} 

h = (z_R - (w_s * ((N / 2) - 1)) - (D*(N/2))) / N; % derive grouser height from wheel diameter and maximum size possible based on size constraint and wheel count
x_g = b; % width of grouser cross sect face {m} (i.e. looking at grouser head-on from the part that touches the ground, this is the width of its face) (equal to wheel width)
N_g = (N_Grouser*l) / ((2*pi)*(D/2)); %number of grousers in contact with the ground (based on arc length of the wheel in contact with the ground and the count of grousers)

% round down if decimal is not at least 0.8 or more
% prety arbitrary but conservative (want 80% contact to say the grouser is in contact with the soil)
if (N_g - floor(N_g)) >= 0.8
    N_g = round(N_g); %round to nearest whole number
else
    N_g = floor(N_g); %round to floor
end

alpha = acosd(1 - ((2*z)/D)); % angle of attack of wheel in soil {deg}

l_o = z*(tand(rad2deg(pi / 4)) - phi/2)^2;

omega = v / ((1 - s)*(D/2)); % omega = (v / (D/2)) * 1.05 % solve for omega from wheel slip ratio

%% Constraint Check

% Length Constraint
length = (D*(N/2) + 2*h*(N/2) + w_s * ((N / 2) - 1));

if (length <= z_R)
    fprintf("PASS: Wheel Diameter, Wheel Count, and Grouser Height are within allowances. %fm is <= z_R = %fm\n", length, z_R);
else
    fprintf("FAIL: Wheel Diameter, Wheel Count, and Grouser Height FAIL the contraint check. One of these is too large. %fm is > z_R = %fm\n", length, z_R);
end

% Weight constraint
W_max_per_wheel = A*(c*N_c + gamma*z*N_q + (1/2)*gamma*b*N_gamma); % max weight on soil per wheel, N
W_max = N*W_max_per_wheel;

if (W < W_max)
    fprintf("PASS: Rover weight and wheel thickness are within allowances. Weight: %f < W_max: %f\n", W, W_max);
else
    fprintf("FAIL: Rover weight and wheel thickness are NOT within allowances. Weight may be too high, or thickness is too small. Weight: %f >= W_max: %f\n", W, W_max);
end

% Grouser cross sectional y_g height or "thickness" constraint 
distance_between_grousers = (2*pi*(D/2))/N_Grouser; % spacing between grousers {m}

if (y_g < distance_between_grousers)
    fprintf("PASS: Grouser y_g cross sectional face height is within allowances. y_g = %fm is < %fm\n", y_g, distance_between_grousers);
else
    fprintf("FAIL: Grouser y_g cross sectional face height is NOT within allowances. y_g = %fm is >= %fm\n", y_g, distance_between_grousers);
end

% Grousers in contact constraint
if (N_g >= 1)
    fprintf("PASS: Number of Grousers in contact with soil: %d is >= 1\n", N_g);
else
    fprintf("FAIL: Number of Grousers in contact with soil: %d is < 1\n", N_g);
end

% **Ignore bulldozing resistance constraint** Based on apostolopoulos_dimitrios_2001_dissertation
ratio = z/D;
if (ratio < 0.06)
    fprintf("PASS: Can ignore Bulldozing Resistance: %f is < 0.06\n", ratio);
else
    fprintf("FAIL: CANNOT ignore Bulldozing Resistance: %f is >= 0.06\n", ratio);
end

%% Drawpull Calculations

R_c = compressionResistance(b, k, n, z);
fprintf("Compression Resistance = " + R_c + " N\n")

R_r = rollingResistance(b, k, W, D);
fprintf("Rolling Resistance = " + R_r + " N\n")

R_g = gravitationalResistance(W, theta);
fprintf("Gravitational Resistance = " + R_g + " N\n")

R_b = bulldozeResistance(alpha, phi, b, z, c, K_bar_c, gamma, K_bar_gamma, l_o, D);
fprintf("Bullzode Resistance = " + R_b + " N\n")

H_smooth = smoothTractiveForce(A, c, W_w, phi, K_bar, s, l, N);
fprintf("Smooth Tractive Force = " + H_smooth + " N\n")

H_grouser = grouserTractiveForce(b, l, c, h, N_g, W_w, phi, K_bar, s, N);
fprintf("Grouser Tractive Force = " + H_grouser + " N\n")

DP_smooth = drawbarPullSmooth(H_smooth, R_c, R_b, R_g, R_r);
fprintf("Drawbar Pull Smooth = " + DP_smooth + " N\n")

DP_grouser = drawbarPullGrouser(H_grouser, R_c, R_b, R_g, R_r);
fprintf("Drawbar Pull Grouser = " + DP_grouser + " N\n")

T_smooth = treqSmooth(H_smooth, N, D);
fprintf("Torque Smooth = " + T_smooth + " N-m\n")

T_grouser = treqGrouser(H_smooth, N, D);
fprintf("Torque Grouser = " + T_grouser + " N-m\n")


%% Grouser Structural Failure

% F1 = H_grouser/N/N_g; % thrust onto grouser, N
% dist = h/2; % equivalent force location for rectangular load, m
% M1 = F1*dist; % moment induced by the load, N-m
% Ix = 1/12*x_g*y_g^3; % moi on grouser in x, m^4
% c_g1 = y_g/2; % height from neutral axis
% sigma_b1 = M1*c_g1/Ix; % bending stress due to pure bending
% 
% F2 = 9.92/N_g; % Side load due to steering (calculated later), N
% M2 = F2*dist; % moment induced by the load, N-m
% Iy = 1/12*y_g*x_g^3; % moi of grouser in x, m^4
% c_g2 = x_g / 2; % height from neutral axis to load {m}
% sigma_b2 = M2*c_g2 / Iy; % bending struess due to side loads
% 
% V1 = (H_grouser/N)/N_g; % thrust load, N
% tau1 = V1 / (y_g*x_g); % force per grouser / cross sect area
% 
% V2 = 9.92/N_g; % side load, N
% tau2 = V2 / (y_g*x_g); % force per grouser / cross sect area
% 
% tauTotal = sqrt(tau1^2 + tau2^2


%% Drawbar Pull Functions

function R_c = compressionResistance(b,k,n,z)
  R_c_one = (b*k/(n+1))*(z^(n+1)) ;

  % compression resistance for entire rover (2 sides to a rover)
  R_c = R_c_one * 2 ;

end

function R_r = rollingResistance(b,k,W,D)
   R_r = (0.580 / ((b*k)^.5)) * (W^(3/2) / D^(3/4));
end

function R_g = gravitationalResistance(W,theta)
   R_g = W*sind(theta);
end

function R_b = bulldozeResistance(alpha, phi, b, z, c, K_bar_c, gamma, K_bar_gamma, l_o, D)
  % Compression resistance per leading wheel {N}
  R_b_leading = ((b * sind(alpha + phi)) / (2 * sind(alpha) * cosd(phi))) * (2 * z * c * K_bar_c + gamma * z^2 * K_bar_gamma) + (l_o^3 * gamma / 3) * ((pi/2) - deg2rad(phi)) + c * l_o^2 * (1 + tand((rad2deg(pi/4) + phi/2)));

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
    H_grouser_per_wheel = (b*l*c*(1+2*h/b)*N_g + W_w*tand(phi)*(1 + 0.64*(h/b)*atan(b/h))) * (1 - (K_bar / (s*l))*(1 - exp((-s*l)/K_bar)));
    H_grouser = N * H_grouser_per_wheel;
end

function DP = drawbarPullSmooth(H_smooth, R_c, R_b, R_g, R_r)
    DP = H_smooth - (R_c + R_b + R_g + R_r);
end

function DP = drawbarPullGrouser(H_grouser, R_c, R_b, R_g, R_r)
    DP = H_grouser - (R_c + R_b + R_g + R_r);
end

function TReq = treqSmooth(H_smooth, N, D)
    TReq = (H_smooth/N)*(D/2);
end

function TReq = treqGrouser(H_grouser, N, D)
    TReq = (H_grouser/N)*(D/2);
end