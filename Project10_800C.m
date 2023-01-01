% Solid oxide Fuel cell modelling
%Rxn @ anode H2 + O-- ---> H2O+2e-
%    @ Cathode 0.5O2 + 2e- ---> O--
% Electrodes are made of Nickel, Electrolyte Yittria Stabilised Zirconia
% Assume electrodes are of same material and at same temperature
% Assume resistance from tri-junction to terminal is negligible

%%Cleanup
clear all
clc
close all

addpath 'D:\Python\matlab\toolbox'

%% Inputs
T = 800+273.15; %K
P = 3e5; %pa 
R = 8.314; %J/mol-K
F = 96485; %C/mol
Tc_H2 = 33.2; %K
Tc_H2O = 647.14; %K
Pc_H2 = 12.83*1.01325; %atm
Pc_H2O = 217.12*1.01325; %atm
M_H2 = 2.016; %g/mol
M_H2O = 18.01528; %g/mol
Pc_N2 = 33.953;
Pc_O2 = 50.5;
Tc_N2 = 126.2;
Tc_O2 = 154.6;
M_N2 = 28.0134;
M_O2 = 31.999;
K0 = 15;
T0 = 1000+273.15;
i_eqc0 = 1000;
K = K0*exp(-(96.5e3/R)*((1/T)-(1/T0)));%15; %S/m
DH2H2O = (((3.64e-4 * ((T/sqrt(Tc_H2*Tc_H2O))^2.334) * (Pc_H2*Pc_H2O)^(1/3) * (Tc_H2*Tc_H2O)^(5/12)*((1/M_H2)+(1/M_H2O))^(1/2))/1)*1e-4)/(P*1e-5); %(3.8378e-3)/3; %m2/s
DO2N2 = (((2.745e-4 * ((T/sqrt(Tc_N2*Tc_O2))^1.823) * (Pc_N2*Pc_O2)^(1/3) * (Tc_N2*Tc_O2)^(5/12)*((1/M_O2)+(1/M_N2))^(1/2))/1)*1e-4)/(P*1e-5); %(2.9417e-4)/3; %m2/s
L_GDLa = 5e-3; %m
L_GDLc = 5e-3; %m
L_YSZ = 50e-6; %m
i_eqc = i_eqc0*exp(-(100e3/R)*((1/T)-(1/T0)));%1000; %A/m2
i_eqa = i_eqc * 100;
Velocity = 1;
Cellnum = 300;
L = 0.5;
w = 10e-3;
h = 5e-3;
A = w*h;
Voltage_Lower = 0.0;
Voltage_Upper = 1.2;
Voltage_Step = 0.01;
%Acell = 0.005/Cellnum;
Acell = L*w/Cellnum;
%% Defining gas
gas = IdealGasMix('GRI30.xml');
iO2 = speciesIndex(gas,'O2');
iN2 = speciesIndex(gas,'N2');
iH2 = speciesIndex(gas,'H2');
iH2O = speciesIndex(gas,'H2O');

nsp = nSpecies(gas);
xanode = zeros(nsp,1);
xcathode = zeros(nsp,1);
%{
xanode = zeros(nsp,1);
xanode(iH2O) = 3/100;
xanode(iH2) = 1-xanode(iH2O);

xcathode = zeros(nsp,1);
xcathode(iO2) = 0.21;%1;
xcathode(iN2) = 0.79;%3.76;
xcathode = xcathode./sum(xcathode);
NH2 = 0.97*(P/(R*T))*A*Velocity;
NH2O = 0.03*(P/(R*T))*A*Velocity;
NO2 = NH2;
NN2 = 3.76*NO2;
%}
NH20 = 0.97*(P/(R*T))*A*Velocity;
NH2O0 = 0.03*(P/(R*T))*A*Velocity;
NO20 = NH20;
NN20 = 3.76*NO20;
Current = 0;
J = 1;
CurrentSum = 0;
for Voltage = Voltage_Lower:Voltage_Step:Voltage_Upper
    Voltage
NH2 = 0.97*(P/(R*T))*A*Velocity;
NH2O = 0.03*(P/(R*T))*A*Velocity;
NO2 = NH2;
NN2 = 3.76*NO2;
xanode(iH2O) = NH2O/(NH2+NH2O);
xanode(iH2) = NH2/(NH2+NH2O);
xcathode(iO2) = NO2/(NO2+NN2);
xcathode(iN2) = NN2/(NO2+NN2);
trigger = 0;
for C = 1:1:Cellnum
    B=1;
    error = 0;
    errorprev = 1;
    errorprev2 = 2;
    i(C) = 5*10^6;
    step = 2.5*10^6;
   
    if C > 1
    i(C) = abs(Current(C-1,J)/Acell);
    step = 0.1*i(C);
    end
    
    if C > 1

set(gas,'T',T,'P',P,'X',xanode)
hanodepre = enthalpy_mole(gas)/1000;
set(gas,'T',T,'P',P,'X',xcathode)
hcathodepre = enthalpy_mole(gas)/1000;
EnthalpyPre = hcathodepre*(NN2 + NO2)/Acell + hanodepre*(NH2+NH2O)/Acell;
NH2 = NH2 - (Current(C-1,J))/(2*F);
NH2O = NH2O + (Current(C-1,J))/(2*F);
NO2 = NO2 - (Current(C-1,J))/(4*F);
NN2 = NN2;

xanode(iH2O) = NH2O/(NH2+NH2O);
xanode(iH2) = NH2/(NH2+NH2O);
xcathode(iO2) = NO2/(NO2+NN2);
xcathode(iN2) = NN2/(NO2+NN2);
%{
cathodedensity = (NO2 + NN2)/(w*h*Velocity);
anodedensity = (NH2+ NH2O)/(Velocity*w*h);
Pcathode = T*(cathodedensity*R);
Panode = T*(anodedensity*R);
%}

set(gas,'T',T,'P',P,'X',xanode)
hanodepost = enthalpy_mole(gas)/1000;
set(gas,'T',T,'P',P,'X',xcathode)
hcathodepost = enthalpy_mole(gas)/1000;

HeatTransferDensity(C-1,J) = ((hcathodepost)*(NN2 + NO2))/Acell + ((hanodepost)*(NH2 + NH2O))/Acell - EnthalpyPre + power(C-1,J)/Acell;
HeatTransfer(C-1,J) = HeatTransferDensity(C-1,J)*Acell;




    end
xH2O(C,J) = xanode(iH2O);
xH2(C,J) = xanode(iH2);
xO2(C,J) = xcathode(iO2);
xN2(C,J) = xcathode(iN2);
Ncathode(C,J) = NO2+NN2;
Nanode(C,J) = NH2+NH2O;
    Error = 1;            
    while Error > 1e-9
%% Equilibrium state calc
%{
if trigger == 1
    i(C) = 0;
    power(C,J) = 0;
    Current(C,J) = 0;
    Error = 1e-8;

end
%}
set(gas,'T',T,'P',P,'X',xanode)
mu_anode = chemPotentials(gas);
mu_h2eqga = mu_anode(iH2)/1000; %J/mol
mu_h2oeqga = mu_anode(iH2O)/1000; %J/mol

set(gas,'T',T,'P',P,'X',xcathode)
mu_cathode = chemPotentials(gas);
mu_o2eqgc = mu_cathode(iO2)/1000; %J/mol

%% Assume
muec_eeqa = 0;

%% Calculate
c = P/(R*T);
muec_oeqysza = mu_h2oeqga + (2*muec_eeqa) - mu_h2eqga;
muec_oeqyszc = muec_oeqysza;
muec_eeqc = (0.5*muec_oeqyszc) - (0.25*mu_o2eqgc);

R_eqa = i_eqa/(2*F);
R_eqc = i_eqc/(2*F);

%ze = -1 zo-- = -2

pot_diff_eq(C,J) = (1/F)*(muec_eeqa - muec_eeqc);
%pot_diff_eq_check = (1/(2*F))*(mu_h2eqga + (0.5*mu_o2eqgc) - mu_h2oeqga)
t = 1;
check = 1;
%i = 0:1e3:50e3;

    v = i(C)/(2*F);
    J_ea = (2*v);
    J_h2a = v;
    J_h2oa = v;
    
    delmuec_ea = 0; % J_ea*L_a/kp_e
    x_h2a(C) = xanode(iH2) - ((J_h2a*L_GDLa)/(c*DH2H2O));
    x_h2oa(C) = xanode(iH2O) + ((J_h2oa*L_GDLa)/(c*DH2H2O));
    
    delmu_h2gdla(C) = -(R*T)*log(1-((J_h2a*L_GDLa)/(xanode(iH2)*c*DH2H2O)));
    delmu_h2ogdla(C) = (R*T)*log(1+((J_h2oa*L_GDLa)/(xanode(iH2O)*c*DH2H2O)));
    
    muec_ea(C) = muec_eeqa + delmuec_ea;
    mu_h2a(C) = mu_h2eqga - delmu_h2gdla(C);
    mu_h2oa(C) = mu_h2oeqga + delmu_h2ogdla(C);
    
    muec_oa(C) = muec_oeqysza + delmu_h2gdla(C) + R*T*log((v/R_eqa)+exp((delmu_h2ogdla(C)+(2*delmuec_ea))/(R*T)));
    
    J_o  = v;
    kp_yszo = K/(4*(F^2));
    delmuec_oysz(C) = J_o*L_YSZ/kp_yszo;
    
    muec_oc(C) = muec_oa(C) + delmuec_oysz(C);
    
    J_o2c = v/2;
    x_o2c(C) = 1- (1-xcathode(iO2))*exp((J_o2c*L_GDLc)/(c*DO2N2));
    
    delmu_o2gdlc(C) = -(R*T)*log(x_o2c(C)/xcathode(iO2));
    
    mu_o2c(C) = mu_o2eqgc - delmu_o2gdlc(C);
    
    muec_ec(C) = muec_eeqc + (0.25*delmu_o2gdlc(C)) + 0.5*R*T*log((v/R_eqc)+exp((muec_oc(C)-muec_oeqyszc)/(R*T)));
    J_ec = 2*v;
    
    delmuec_ec = 0; %J_ec*L_c/kp_e
    
    muec_etc = muec_ec(C) + delmuec_ec;
    muec_eta = muec_eeqa;
    pot_difft(C) = -(1/F)*(muec_etc-muec_eta);
    power(C,J) = pot_difft(C)*i(C)*Acell;
    Current(C,J) = i(C)*Acell;
    if B > 2
        errorprev2 = errorprev;
    end
    if B > 1
    errorprev = error;
    end
    Error = abs(pot_difft(C)-Voltage);

    error = pot_difft(C)-Voltage;
    if error > 0 && x_o2c(C) > 0
        i(C) = i(C) + step;
        step = step/2;
    elseif error < 0 && x_o2c(C) > 0
        i(C) = i(C) - step;
        step = step/2;
    end
    if x_o2c(C) < 0 || x_h2oa(C) > 1 || x_h2a(C) < 0
        i(C) = i(C) - step;
        step = step/2;
    end
    power(C,J) = pot_difft(C)*i(C)*Acell;
    
    B = B+1;
    if B > 25 && errorprev == error && errorprev == errorprev2 && Error > 1e-4 || B > 100
        i(C) = 2*F*NH2/Acell;
        Error = 1e-10;
        
        
    end
    end
if C > 1
    CurrentSum(C,J) = CurrentSum(C-1,J) + Current(C,J);
else
    CurrentSum(C,J) = Current(C,J);
end

end
Volt(J) = Voltage;
J = J+1;
end
    
    
%% 
figure
hold on
totalpower = sum(power,1);
totalcurrent = sum(Current,1);
totalheatTransfer = -sum(HeatTransfer,1);
FuelUtil = 0.97 - xH2(end,:);
plot(totalcurrent,totalpower/100)
plot(totalcurrent,Volt)
plot(totalcurrent,FuelUtil)
plot(totalcurrent,totalheatTransfer/100)
legend('Power Output/100 (W)','Electrical Potential (V)','Fuel Utilization','Heat Transfer/100 (W)')
xlabel('Cell Current (A)')
title('800 C, 3 bar, 1 m/s')
xlim([0 200])
ylim([0 1.75])
hold off
plotfixer



figure
hold on
plot(totalcurrent,Volt)
LHV_hydrogen = 119.96e6; %J/kg
LHV_hydrogen_out = xH2(end,:).*Nanode(end,:)*LHV_hydrogen*0.002;
LHV_hydrogen_in = xH2(1,:).*Nanode(1,:)*LHV_hydrogen*0.002;
plot(totalcurrent,LHV_hydrogen_out/200)
plot(totalcurrent,totalpower./LHV_hydrogen_in)
plot(totalcurrent,totalpower./(LHV_hydrogen_in-LHV_hydrogen_out))
legend('Electrical Potential (V)','H2 Heating Value/200 (W)','Efficiency (LHV, H2in)','Efficiency (LHV, H2used)')
xlabel('Cell Current (A)')
title('800 C, 3 bar, 1 m/s')
xlim([0 200])
ylim([0 2])
hold off
plotfixer

figure
hold on
Distance = L*(1:1:Cellnum)/Cellnum;
plot([0 Distance],[0;(Current(:,55)/Acell)]/1000)
plot([0 Distance],[0;(power(:,55)/Acell)]/1000)
plot([0 Distance(1:end-1)],[0;(-HeatTransfer(:,55)/Acell)]/1000)
legend('Current Density (kA/m^{2})','Elec. Power Density (kW/m^{2})','Heat Flux (kW/m^{2})')
title('800 C, 3 bar, 0.54 V')
xlabel('Distance Along Channel (m)')
xlim([0 0.5])
ylim([0 30])
hold off
plotfixer

figure
hold on
Distance = L*(1:1:Cellnum)/Cellnum;
plot(Distance,xH2(:,55))
plot(Distance,xO2(:,55))
plot(Distance,xH2O(:,55))
plot(Distance,pot_diff_eq(:,55))
plot([0 0.5],[0.54 0.54],'--k')
legend('xH2 (anode)','xO2 (cathode)','xH2O (anode)','\Delta \phi equil.(V)','\Delta \phi cell (V)')
title('800 C, 3 bar, 0.54 V')
xlabel('Distance Along Channel (m)')
xlim([0 0.5])
ylim([0 1.5])
hold off
plotfixer


%{
figure
hold on
Distance = L*(1:1:Cellnum)/Cellnum;
plot([0 Distance],[0;(Current(:)/Acell)]/1000)
plot([0 Distance],[0;(power(:)/Acell)]/1000)
plot([0 Distance(1:end-1)],[0;(-HeatTransfer(:)/Acell)]/1000)
legend('Current Density (kA/m^{2})','Elec. Power Density (kW/m^{2})','Heat Flux (kW/m^{2})')
title('1000 C, 3 bar, 0.50 V, Velocity = 0.5 m/s')
xlabel('Distance Along Channel (m)')
xlim([0 0.5])
ylim([0 60])
hold off
plotfixer

figure
hold on
Distance = L*(1:1:Cellnum)/Cellnum;
plot(Distance,xH2(:))
plot(Distance,xO2(:))
plot(Distance,xH2O(:))
plot(Distance,pot_diff_eq(:))
plot([0 0.5],[0.50 0.50],'--k')
legend('xH2 (anode)','xO2 (cathode)','xH2O (anode)','\Delta \phi equil.(V)','\Delta \phi cell (V)')
title('1000 C, 3 bar, 0.50 V, Velocity = 0.5 m/s')
xlabel('Distance Along Channel (m)')
xlim([0 0.5])
ylim([0 1.5])
hold off
plotfixer
%}












%{
figure(1)
plot(i/1000,pot_difft,'b')
hold on
plot(i/1000,power/max(power),'g')
plot(i/1000,x_h2a,'r')
plot(i/1000,x_h2oa,'m')
plot(i/1000,x_o2c,'k')
plot([0 60],[0.97 0.97],'--r')
plot([0 60],[0.21 0.21],'--k')
plot([0 60],[0.03 0.03],'--m')
xlabel('Current density [kA/m2]')
ylabel('')
title('YSZ: 1000C 1bar')
legend('Electric potential(V)','Power/Max Power','Anode Hydrogen','Anode Water','Cathode oxygen')
str1 = {'Max power is '};
str2 = num2str(max(power)/1000);
str3 = {' kW'};
str = strcat(str1,str2,str3);
text(10,1.05,str)
%xlim([0 60])
ylim([0 1.2])
    
Gas_diffusion_loss = (delmu_h2gdla) + (delmu_h2ogdla) + (0.5*delmu_o2gdlc);
Gas_diffusion_loss = Gas_diffusion_loss/(96485); %What is mol-rxn
Ohmic_loss = (delmuec_oysz)/(96485);
Activation_loss_anode = (-mu_h2oa - (2*muec_ea) + mu_h2a + muec_oa)/(96485);
Activation_loss_cathode = (-muec_oc + (0.5*mu_o2c) + (2*muec_ec))/(96485);
figure(2)
plot(i/1000,Gas_diffusion_loss,'b')
hold on
plot(i/1000,Ohmic_loss,'r')
plot(i/1000,Activation_loss_cathode,'k')
plot(i/1000,Activation_loss_anode,'m')
legend('GDL Loss','Ohmic Loss','Cathode loss','Anode loss')
xlabel('Current density [kA/m2]')
ylabel('Electrochemical potential[eV/rxn]')
title('YSZ: 1000C 1bar')
xlim([0 60])
%ylim([0 0.7])

overall_affinity = mu_h2eqga + (0.5*mu_o2eqgc) - mu_h2oeqga;
ov_affinity = overall_affinity*ones(1,length(i))/96485;
figure(3)
plot(i/1000,ov_affinity,'b')
hold on
% ov_GDL = ov_affinity - Gas_diffusion_loss;
% plot(i/1000,ov_GDL)
% ov_ohmic = ov_GDL - Ohmic_loss;
% plot(i/1000,ov_ohmic)
% ov_anode = ov_ohmic - Activation_loss_anode;
% plot(i/1000,ov_anode)
% ov_cathode = ov_anode - Activation_loss_cathode;
% plot(i/1000,ov_cathode)

ov_anode = ov_affinity - Activation_loss_anode;
plot(i/1000,ov_anode,'k')
ov_cathode = ov_anode - Activation_loss_cathode;
plot(i/1000,ov_cathode,'g')
ov_ohmic = ov_cathode - Ohmic_loss;
plot(i/1000,ov_ohmic,'m')
ov_GDL = ov_ohmic - Gas_diffusion_loss;
plot(i/1000,ov_GDL,'r')
xlabel('Current density [kA/m2]')
ylabel('Electrochemical potential[eV/rxn]')
title('Impact of losses')
legend('Overall Affinity','- Anode loss','- Cathode loss','- Ohmic loss','- Gas diffusion loss')
xlim([0 60])
%}