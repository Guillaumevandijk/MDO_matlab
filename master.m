clear all;
close all;

params = refParams();

% Wing planform geometry 
%                x    y     z   chord(m)    twist angle (deg) 
AC.Wing.Geom = [0     0     0     6.05         0;
                2.26  4.73  0.698     4.16         0;
                5.76  14   1.184     1.42         0];

% Wing incidence angle (degree)
AC.Wing.inc  = 0;
kinkLoc = params.kinkLoc;
Au1 = [0.2171    0.3450    0.2975    0.2685    0.2893];
Al1 = [-0.1299   -0.2388   -0.1635   -0.0476    0.0797];

Au3 = [0.2171     0.4    0.2975    0.2685    0.2893];
Al3 = [-0.1299   -0.15   -0.1635   -0.0476    0.0797];

Au2 = (Au1*(1-kinkLoc)+Au3*kinkLoc);
Al2 = (Al1*(1-kinkLoc)+Al3*kinkLoc);
            
% Airfoil coefficients input matrix
%                    | ->     upper curve coeff.                <-|   | ->       lower curve coeff.       <-| 
AC.Wing.Airfoils   =    [Au1  Al1;
                        Au2  Al2
                        Au3  Al3];
                  

AC.Wing.eta = [0;kinkLoc;1];  % Spanwise location of the airfoil sections

% Viscous vs inviscid
AC.Visc  = 0;              % 0 for inviscid and 1 for viscous analysis
AC.Aero.MaxIterIndex = 150;    %Maximum number of Iteration for the
                                %convergence of viscous calculation
                 

V_maxCruise = params.V_maxCruise;
meanChord = params.meanChord;
rho = params.density_cruise;
alt = params.alt_cruise;
viscosity = params.viscosity_cruise;
T_cruise = params.T_cruise;
Re = rho* meanChord*V_maxCruise/viscosity;
AoA = params.AoA;
a = sqrt(1.4*287*T_cruise);

                                
% Flight Condition
AC.Aero.V     = V_maxCruise;            % flight speed (m/s)
AC.Aero.rho   = rho;         % air density  (kg/m3)
AC.Aero.alt   = alt;             % flight altitude (m)
AC.Aero.Re    = Re;        % reynolds number (bqased on mean aerodynamic chord)
AC.Aero.M     = V_maxCruise/a;           % flight Mach number 
% AC.Aero.CL    = 0.4;          % lift coefficient - comment this line to run the code for given alpha%
AC.Aero.Alpha = params.AoA;             % angle of attack -  comment this line to run the code for given cl 


%% 
tic

Res = Q3D_solver(AC);
disp(Res.Wing.cl);
toc

xpoints = linspace(0,1,50)';
[Xtu1,Xtl1,C1] = D_airfoil2(Au1,Al1,xpoints);
[Xtu2,Xtl2,C2] = D_airfoil2(Au2,Al2,xpoints);
[Xtu3,Xtl3,C3] = D_airfoil2(Au3,Al3,xpoints);

% A=importdata('e553.dat');
% upper= A(1:35,:);
% lower = A(36:73,:);
% lower(1,:) = [0 0]; 
% upper = [0,0;upper];
% upper_Out = interp1(upper(:,1),upper(:,2),xpoints);
% lower_Out = interp1(lower(:,1),lower(:,2),xpoints);


% hold on
% plot(Xtu1(:,1),Xtu1(:,2),'b');    %plot upper surface coords
% plot(Xtl1(:,1),Xtl1(:,2),'b');    %plot lower surface coords
% plot(Xtu2(:,1),Xtu2(:,2),'r');    %plot upper surface coords
% plot(Xtl2(:,1),Xtl2(:,2),'r');    %plot lower surface coords
% plot(Xtu3(:,1),Xtu3(:,2),'g');    %plot upper surface coords
% plot(Xtl3(:,1),Xtl3(:,2),'g');    %plot lower surface coords
% %plot(xpoints,lower_Out,'r');
% %plot(xpoints,upper_Out,'r');
% axis([0,1,-.5,.5]);
% hold on;