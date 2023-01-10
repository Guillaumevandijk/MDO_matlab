clear all;
close all;

params = refParams();
global loc_kink;
loc_kink = params.loc_kink;


%Defign the root and tip airfoil

Au1 = [0.2171    0.3450    0.2975    0.2685    0.2893];
Al1 = [-0.1299   -0.2388   -0.1635   -0.0476    0.0797];

Au3 = [0.2171     0.4    0.2975    0.2685    0.2893];
Al3 = [-0.1299   -0.15   -0.1635   -0.0476    0.0797];

%Calculate the airfoil at the kink by linear interpolation between the root
%and tip chord.

Au2 = (Au1*(1-loc_kink)+Au3*loc_kink);
Al2 = (Al1*(1-loc_kink)+Al3*loc_kink);

%put these coefficients in the design vector

designVector.coeffs = [Au1  Al1;
                        Au2  Al2
                        Au3  Al3];

designVector.chord1 = 6.05;
designVector.chord2 = 4.16;
designVector.chord3 = 1.42;
designVector.span = 28.076;


span = designVector.span;




disp(constraints(designVector));
res = loads(designVector);







%% 
% tic
% 
% Res = Q3D_solver(AC);
% disp(Res.Wing.cl);
% toc


