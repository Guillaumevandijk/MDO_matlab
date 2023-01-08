clear all;
close all;

hfCO2 =     -94.01;     %kcal/mol
hfH20 =     -57.78;     %kcal/mol
hfC12H24 =  -92.2;      %kcal/mol
hc_gram    =     -10.278;    %kcal/gram

molarWeight = 168;      %gram/mol
hc_mol = hc_gram*168;

disp("hc per mol =");
disp(hc_mol);


hc_calc = 12*hfCO2+12*hfH20-hfC12H24;
disp("calculated hc per mol=")
disp(hc_calc);



