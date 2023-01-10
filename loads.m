function [loadVector] = loads(x)

coeffs = x.coeffs;
span = x.span;
chord1 = x.chord1;
chord2 = x.chord2;
chord3 = x.chord3;

coeffsNr = 5;

Au1 = coeffs(1,1:coeffsNr);
Al1 = coeffs(1,coeffsNr+1:coeffsNr*2);
Au2 = coeffs(2,1:coeffsNr);
Al2 = coeffs(2,coeffsNr+1:coeffsNr*2);
Au3 = coeffs(3,1:coeffsNr);
Al3 = coeffs(3,coeffsNr+1:coeffsNr*2);
params = refParams();

global loc_kink;

% Wing planform geometry 
%                x    y     z   chord(m)    twist angle (deg) 
AC.Wing.Geom = [0     0     0     chord1         0;
                2.26  4.73  0.698     chord2         0;
                5.76  span/2   1.184     chord3         0];

% Wing incidence angle (degree)
AC.Wing.inc  = 0;
            
% Airfoil coefficients input matrix
%                    | ->     upper curve coeff.                <-|   | ->       lower curve coeff.       <-| 
AC.Wing.Airfoils   =    [Au1 Al1;
                        Au2 Al2;
                        Au3 Al3];
                  

AC.Wing.eta = [0;loc_kink;1];  % Spanwise location of the airfoil sections

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

loadVector = [Res.Wing.Yst Res.Wing.cl]; 



end