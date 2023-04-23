clc; clear; close all;

% Initialize values for global variables
m = 50; % rover mass (max), kg
b = 4/39.37; % width of wheels, m
D = 10/39.37; % wheel diameter, m
N = 6; % number of wheels
w_s = 0.125; % wheel spacing from edge of grouser to edge of grouser {m} 
g = 9.81; % gravity on earth {m/s^2}

% Dynamics
v = 1.5; % rover speed {m/s}*
s = 0.5; % s = 1 - v / (omega * (D/2)) % wheel slip ratio (assume a conservative value of 0.5)

% Grousers
y_g = 0.25/39.37; % height of grouser cross sect face {m} (i.e. looking at grouser head-on from the part that touches the ground, this is the height of its face)
N_Grouser = 60; % total number of grousers

wheel_calcs;

%% Wheel Diameter Optimization
DBP = [];
diam = linspace(3.5, 12, 1000);
for i = diam
    D = i/39.37;
    wheel_calcs;
    DBP = [-DP_grouser*8.8507457676 DBP];
end

figure
plot(diam, DBP)
title("Drawbar Pull vs Diameter")
xlabel("Wheel Diameter (in)")
ylabel("Drawbar Pull (lb)")

%% Wheel Width Optimization
% DBP = [];
% width = linspace(3, 12, 1000);
% for i = width
%     b = i/39.37;
%     wheel_calcs;
%     DBP = [-DP_grouser*8.8507457676 DBP];
% end
% 
% figure
% plot(width, DBP)
% title("Drawbar Pull vs Wdith")
% xlabel("Wheel Width (in)")
% ylabel("Drawbar Pull (lb)")

%% Wheelspan Optimization
% DBP = [];
% wheelspan = linspace(0.01, 0.14, 1000);
% for i = wheelspan
%     w_s = i/39.37;
%     wheel_calcs;
%     DBP = [-DP_grouser*8.8507457676 DBP];
% end
% 
% figure
% plot(wheelspan, DBP)
% title("Drawbar Pull vs Wheelspan")
% xlabel("Distance Between Grousers (in)")
% ylabel("Drawbar Pull (lb)")

%% Grousers Optimization
% DBP = [];
% grousers = linspace(20, 60, 1000);
% for i = grousers
%     N_Grouser = i;
%     wheel_calcs;
%     DBP = [-DP_grouser*8.8507457676 DBP];
% end
% 
% figure
% plot(grousers, DBP)
% title("Drawbar Pull vs Number of Grousers")
% xlabel("Number of Grousers")
% ylabel("Drawbar Pull (lb)")

%% Mass Optimization
% DBP = [];
% mass = linspace(20, 50, 1000);
% for i = mass
%     m = i;
%     wheel_calcs;
%     DBP = [-DP_grouser*8.8507457676 DBP];
% end
% 
% figure
% plot(mass, DBP)
% title("Drawbar Pull vs Mass")
% xlabel("Mass (kg)")
% ylabel("Drawbar Pull (lb)")

%% Slip Ratio Optimization
% DBP = [];
% slip = linspace(-1, 1, 1000);
% for i = slip
%     s = i;
%     wheel_calcs;
%     DBP = [-DP_grouser*8.8507457676 DBP];
% end
% 
% figure
% plot(slip, DBP)
% title("Drawbar Pull vs Slip")
% xlabel("Slip")
% ylabel("Drawbar Pull (lb)")

