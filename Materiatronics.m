%% Definitions
clear all
close all
clc

name='example_file'; % put here the name of your .csv file from ANSYS HFSS (v. 2017 and newer). Do not include ".csv" in the name.
power=-3; % paste here the units along the Y-axis of the plot from the ANSYS HFSS project: Volt(0), miliVolt(-3), mikroVolt(-6), nanoVolt(-9), or pikoVolt(-12), etc.
fpower=9; % paste here the units along the X-axis of the plot from the ANSYS HFSS project: MHz(6), GHz(9), THz(12), etc.
fcentral=7.9e9; % paste here the frequency in Hz at which the modular decomposition will be applied (within your frequency sweep)
P=importdata([name '.csv']);
P=P.data.*10^(power);
f=10^(fpower-power).*P(:,1); E0=1;
for ff=1:length(f)
    delta(ff)=abs(f(ff)-fcentral);
end;
[mindelta,f0] = min(delta);
'The frequency at which decomposition was performed is:'
f_central=f(f0)
c=299792458; m0=4.*pi.*10^(-7); e0=1./(m0.*c^2);
gamma=(pi.*f.^2)/(e0.*c^2); eta=sqrt(m0./e0);

set(0,'DefaultAxesFontName', 'Times New Roman')
set(0,'DefaultAxesFontSize', 14)

%% Importing data
%incidence xy
pzXex_hy=P(:,2)+1j.*P(:,3);    mzXex_hy=P(:,4)+1j.*P(:,5);    pzXex_nhy=P(:,6)+1j.*P(:,7);    mzXex_nhy=P(:,8)+1j.*P(:,9);
pzYex_hy=P(:,10)+1j.*P(:,11);  mzYex_hy=P(:,12)+1j.*P(:,13);  pzYex_nhy=P(:,14)+1j.*P(:,15);  mzYex_nhy=P(:,16)+1j.*P(:,17);
pzXey_hx=P(:,18)+1j.*P(:,19);  mzXey_hx=P(:,20)+1j.*P(:,21);  pzXey_nhx=P(:,22)+1j.*P(:,23);  mzXey_nhx=P(:,24)+1j.*P(:,25);
pzYey_hx=P(:,26)+1j.*P(:,27);  mzYey_hx=P(:,28)+1j.*P(:,29);  pzYey_nhx=P(:,30)+1j.*P(:,31);  mzYey_nhx=P(:,32)+1j.*P(:,33);
%incidence xz
pyXex_hz=P(:,34)+1j.*P(:,35);  myXex_hz=P(:,36)+1j.*P(:,37);  pyXex_nhz=P(:,38)+1j.*P(:,39);  myXex_nhz=P(:,40)+1j.*P(:,41);
pyZex_hz=P(:,42)+1j.*P(:,43);  myZex_hz=P(:,44)+1j.*P(:,45);  pyZex_nhz=P(:,46)+1j.*P(:,47);  myZex_nhz=P(:,48)+1j.*P(:,49);
pyXez_hx=P(:,50)+1j.*P(:,51);  myXez_hx=P(:,52)+1j.*P(:,53);  pyXez_nhx=P(:,54)+1j.*P(:,55);  myXez_nhx=P(:,56)+1j.*P(:,57);
pyZez_hx=P(:,58)+1j.*P(:,59);  myZez_hx=P(:,60)+1j.*P(:,61);  pyZez_nhx=P(:,62)+1j.*P(:,63);  myZez_nhx=P(:,64)+1j.*P(:,65);
%incidence yz
pxYey_hz=P(:,66)+1j.*P(:,67);  mxYey_hz=P(:,68)+1j.*P(:,69);  pxYey_nhz=P(:,70)+1j.*P(:,71);  mxYey_nhz=P(:,72)+1j.*P(:,73);
pxZey_hz=P(:,74)+1j.*P(:,75);  mxZey_hz=P(:,76)+1j.*P(:,77);  pxZey_nhz=P(:,78)+1j.*P(:,79);  mxZey_nhz=P(:,80)+1j.*P(:,81);
pxYez_hy=P(:,82)+1j.*P(:,83);  mxYez_hy=P(:,84)+1j.*P(:,85);  pxYez_nhy=P(:,86)+1j.*P(:,87);  mxYez_nhy=P(:,88)+1j.*P(:,89);
pxZez_hy=P(:,90)+1j.*P(:,91);  mxZez_hy=P(:,92)+1j.*P(:,93);  pxZez_nhy=P(:,94)+1j.*P(:,95);  mxZez_nhy=P(:,96)+1j.*P(:,97);

%% Calculating polarizability tensors

%%%%%%%%%% incidence xy
ee11=(1./(4.*gamma.*E0)).*     ( pzXex_hy + mzXex_hy + pzXex_nhy + mzXex_nhy );
em12=(eta./(4.*gamma.*E0)).*   ( pzXex_hy + mzXex_hy - pzXex_nhy - mzXex_nhy );
ee21=(1./(4.*gamma.*E0)).*     ( pzYex_hy + mzYex_hy + pzYex_nhy + mzYex_nhy );
em22=(eta./(4.*gamma.*E0)).*   ( pzYex_hy + mzYex_hy - pzYex_nhy - mzYex_nhy );
me11=(eta./(4.*gamma.*E0)).*   ( mzYex_hy - pzYex_hy + mzYex_nhy - pzYex_nhy );
mm12=(eta^2./(4.*gamma.*E0)).* ( mzYex_hy - pzYex_hy - mzYex_nhy + pzYex_nhy );
me21=(eta./(4.*gamma.*E0)).*   ( pzXex_hy - mzXex_hy + pzXex_nhy - mzXex_nhy );
mm22=(eta^2./(4.*gamma.*E0)).* ( pzXex_hy - mzXex_hy - pzXex_nhy + mzXex_nhy );

ee12=(1./(4.*gamma.*E0)).*     ( pzXey_hx + mzXey_hx + pzXey_nhx + mzXey_nhx );
em11=(eta./(4.*gamma.*E0)).*   ( pzXey_hx + mzXey_hx - pzXey_nhx - mzXey_nhx );
ee22=(1./(4.*gamma.*E0)).*     ( pzYey_hx + mzYey_hx + pzYey_nhx + mzYey_nhx );
em21=(eta./(4.*gamma.*E0)).*   ( pzYey_hx + mzYey_hx - pzYey_nhx - mzYey_nhx );
me12=(eta./(4.*gamma.*E0)).*   ( mzYey_hx - pzYey_hx + mzYey_nhx - pzYey_nhx );
mm11=(eta^2./(4.*gamma.*E0)).* ( mzYey_hx - pzYey_hx - mzYey_nhx + pzYey_nhx );
me22=(eta./(4.*gamma.*E0)).*   ( pzXey_hx - mzXey_hx + pzXey_nhx - mzXey_nhx );
mm21=(eta^2./(4.*gamma.*E0)).* ( pzXey_hx - mzXey_hx - pzXey_nhx + mzXey_nhx );


%%%%%%%%%% incidence xz
ee11copy=(1./(4.*gamma.*E0)).*     ( pyXex_hz + myXex_hz + pyXex_nhz + myXex_nhz );
em13=(eta./(4.*gamma.*E0)).*       ( pyXex_hz + myXex_hz - pyXex_nhz - myXex_nhz );
ee31=(1./(4.*gamma.*E0)).*         ( pyZex_hz + myZex_hz + pyZex_nhz + myZex_nhz );
em33=(eta./(4.*gamma.*E0)).*       ( pyZex_hz + myZex_hz - pyZex_nhz - myZex_nhz );
me11copy=(eta./(4.*gamma.*E0)).*   ( pyZex_hz - myZex_hz + pyZex_nhz - myZex_nhz );
mm13=(eta^2./(4.*gamma.*E0)).*     ( pyZex_hz - myZex_hz - pyZex_nhz + myZex_nhz );
me31=(eta./(4.*gamma.*E0)).*       ( myXex_hz - pyXex_hz + myXex_nhz - pyXex_nhz );
mm33=(eta^2./(4.*gamma.*E0)).*     ( myXex_hz - pyXex_hz - myXex_nhz + pyXex_nhz );

ee13=(1./(4.*gamma.*E0)).*         ( pyXez_hx + myXez_hx + pyXez_nhx + myXez_nhx );
em11copy=(eta./(4.*gamma.*E0)).*   ( pyXez_hx + myXez_hx - pyXez_nhx - myXez_nhx );
ee33=(1./(4.*gamma.*E0)).*         ( pyZez_hx + myZez_hx + pyZez_nhx + myZez_nhx );
em31=(eta./(4.*gamma.*E0)).*       ( pyZez_hx + myZez_hx - pyZez_nhx - myZez_nhx );
me13=(eta./(4.*gamma.*E0)).*       ( pyZez_hx - myZez_hx + pyZez_nhx - myZez_nhx );
mm11copy=(eta^2./(4.*gamma.*E0)).* ( pyZez_hx - myZez_hx - pyZez_nhx + myZez_nhx );
me33=(eta./(4.*gamma.*E0)).*       ( myXez_hx - pyXez_hx + myXez_nhx - pyXez_nhx );
mm31=(eta^2./(4.*gamma.*E0)).*     ( myXez_hx - pyXez_hx - myXez_nhx + pyXez_nhx );

%%%%%%%%%% incidence yz
ee22copy=(1./(4.*gamma.*E0)).*     ( pxYey_hz + mxYey_hz + pxYey_nhz + mxYey_nhz );
em23=(eta./(4.*gamma.*E0)).*       ( pxYey_hz + mxYey_hz - pxYey_nhz - mxYey_nhz );
ee32=(1./(4.*gamma.*E0)).*         ( pxZey_hz + mxZey_hz + pxZey_nhz + mxZey_nhz );
em33copy=(eta./(4.*gamma.*E0)).*   ( pxZey_hz + mxZey_hz - pxZey_nhz - mxZey_nhz );
me22copy=(eta./(4.*gamma.*E0)).*   ( mxZey_hz - pxZey_hz + mxZey_nhz - pxZey_nhz );
mm23=(eta^2./(4.*gamma.*E0)).*     ( mxZey_hz - pxZey_hz - mxZey_nhz + pxZey_nhz );
me32=(eta./(4.*gamma.*E0)).*       ( pxYey_hz - mxYey_hz + pxYey_nhz - mxYey_nhz );
mm33copy=(eta^2./(4.*gamma.*E0)).* ( pxYey_hz - mxYey_hz - pxYey_nhz + mxYey_nhz );

ee23=(1./(4.*gamma.*E0)).*         ( pxYez_hy + mxYez_hy + pxYez_nhy + mxYez_nhy );
em22copy=(eta./(4.*gamma.*E0)).*   ( pxYez_hy + mxYez_hy - pxYez_nhy - mxYez_nhy );
ee33copy=(1./(4.*gamma.*E0)).*     ( pxZez_hy + mxZez_hy + pxZez_nhy + mxZez_nhy );
em32=(eta./(4.*gamma.*E0)).*       ( pxZez_hy + mxZez_hy - pxZez_nhy - mxZez_nhy );
me23=(eta./(4.*gamma.*E0)).*       ( mxZez_hy - pxZez_hy + mxZez_nhy - pxZez_nhy );
mm22copy=(eta^2./(4.*gamma.*E0)).* ( mxZez_hy - pxZez_hy - mxZez_nhy + pxZez_nhy );
me33copy=(eta./(4.*gamma.*E0)).*   ( pxYez_hy - mxYez_hy + pxYez_nhy - mxYez_nhy );
mm32=(eta^2./(4.*gamma.*E0)).*     ( pxYez_hy - mxYez_hy - pxYez_nhy + mxYez_nhy );

f=f/10^(fpower);  % only for plots

%% Plotting 9 components of alpha_ee tensor

figure8 = figure('Name','9 components of alpha_ee tensor');
set(0,'defaultlinelinewidth',3)
subplot(3,3,1);
hold all
plot(f,eta.*real(ee11),'r')
plot(f,eta.*imag(ee11),'b')
hold off; xlabel(['Frequency [10^{' num2str(fpower) '} Hz]'])
grid on; ylabel('$$\eta_0 \alpha_{\rm ee}^{11} \rm{[m^2 s]}$$','Interpreter','latex')
legend_handle=legend('$$\Re$$','$$\Im$$'); set(legend_handle,'Interpreter','latex');
subplot(3,3,2);
hold all
plot(f,eta.*real(ee12),'r')
plot(f,eta.*imag(ee12),'b')
hold off; xlabel(['Frequency [10^{' num2str(fpower) '} Hz]'])
grid on; ylabel('$$\eta_0 \alpha_{\rm ee}^{12} \rm{[m^2 s]}$$','Interpreter','latex')
legend_handle=legend('$$\Re$$','$$\Im$$'); set(legend_handle,'Interpreter','latex');
subplot(3,3,3);
hold all
plot(f,eta.*real(ee13),'r')
plot(f,eta.*imag(ee13),'b')
hold off; xlabel(['Frequency [10^{' num2str(fpower) '} Hz]'])
grid on; ylabel('$$\eta_0 \alpha_{\rm ee}^{13} \rm{[m^2 s]}$$','Interpreter','latex')
legend_handle=legend('$$\Re$$','$$\Im$$'); set(legend_handle,'Interpreter','latex');
subplot(3,3,4);
hold all
plot(f,eta.*real(ee21),'r')
plot(f,eta.*imag(ee21),'b')
hold off; xlabel(['Frequency [10^{' num2str(fpower) '} Hz]'])
grid on; ylabel('$$\eta_0 \alpha_{\rm ee}^{21} \rm{[m^2 s]}$$','Interpreter','latex')
legend_handle=legend('$$\Re$$','$$\Im$$'); set(legend_handle,'Interpreter','latex');
subplot(3,3,5);
hold all
plot(f,eta.*real(ee22),'r')
plot(f,eta.*imag(ee22),'b')
hold off; xlabel(['Frequency [10^{' num2str(fpower) '} Hz]'])
grid on; ylabel('$$\eta_0 \alpha_{\rm ee}^{22} \rm{[m^2 s]}$$','Interpreter','latex')
legend_handle=legend('$$\Re$$','$$\Im$$'); set(legend_handle,'Interpreter','latex');
subplot(3,3,6);
hold all
plot(f,eta.*real(ee23),'r')
plot(f,eta.*imag(ee23),'b')
hold off; xlabel(['Frequency [10^{' num2str(fpower) '} Hz]'])
grid on; ylabel('$$\eta_0 \alpha_{\rm ee}^{23} \rm{[m^2 s]}$$','Interpreter','latex')
legend_handle=legend('$$\Re$$','$$\Im$$'); set(legend_handle,'Interpreter','latex');
subplot(3,3,7);
hold all
plot(f,eta.*real(ee31),'r')
plot(f,eta.*imag(ee31),'b')
hold off; xlabel(['Frequency [10^{' num2str(fpower) '} Hz]'])
grid on; ylabel('$$\eta_0 \alpha_{\rm ee}^{31} \rm{[m^2 s]}$$','Interpreter','latex')
legend_handle=legend('$$\Re$$','$$\Im$$'); set(legend_handle,'Interpreter','latex');
subplot(3,3,8);
hold all
plot(f,eta.*real(ee32),'r')
plot(f,eta.*imag(ee32),'b')
hold off; xlabel(['Frequency [10^{' num2str(fpower) '} Hz]'])
grid on; ylabel('$$\eta_0 \alpha_{\rm ee}^{32} \rm{[m^2 s]}$$','Interpreter','latex')
legend_handle=legend('$$\Re$$','$$\Im$$'); set(legend_handle,'Interpreter','latex');
subplot(3,3,9);
hold all
plot(f,eta.*real(ee33),'r')
plot(f,eta.*imag(ee33),'b')
hold off; xlabel(['Frequency [10^{' num2str(fpower) '} Hz]'])
grid on; ylabel('$$\eta_0 \alpha_{\rm ee}^{33} \rm{[m^2 s]}$$','Interpreter','latex')
legend_handle=legend('$$\Re$$','$$\Im$$'); set(legend_handle,'Interpreter','latex');

%% Plotting 9 components of alpha_em tensor

figure9 = figure('Name','9 components of alpha_em tensor');
set(0,'defaultlinelinewidth',3)
subplot(3,3,1);
hold all
plot(f, real(em11),'r')
plot(f, imag(em11),'b')
hold off; xlabel(['Frequency [10^{' num2str(fpower) '} Hz]'])
grid on; ylabel('$$ \alpha_{\rm em}^{11} \rm{[m^2 s]}$$','Interpreter','latex')
legend_handle=legend('$$\Re$$','$$\Im$$'); set(legend_handle,'Interpreter','latex');
subplot(3,3,2);
hold all
plot(f, real(em12),'r')
plot(f, imag(em12),'b')
hold off; xlabel(['Frequency [10^{' num2str(fpower) '} Hz]'])
grid on; ylabel('$$ \alpha_{\rm em}^{12} \rm{[m^2 s]}$$','Interpreter','latex')
legend_handle=legend('$$\Re$$','$$\Im$$'); set(legend_handle,'Interpreter','latex');
subplot(3,3,3);
hold all
plot(f, real(em13),'r')
plot(f, imag(em13),'b')
hold off; xlabel(['Frequency [10^{' num2str(fpower) '} Hz]'])
grid on; ylabel('$$ \alpha_{\rm em}^{13} \rm{[m^2 s]}$$','Interpreter','latex')
legend_handle=legend('$$\Re$$','$$\Im$$'); set(legend_handle,'Interpreter','latex');
subplot(3,3,4);
hold all
plot(f, real(em21),'r')
plot(f, imag(em21),'b')
hold off; xlabel(['Frequency [10^{' num2str(fpower) '} Hz]'])
grid on; ylabel('$$ \alpha_{\rm em}^{21} \rm{[m^2 s]}$$','Interpreter','latex')
legend_handle=legend('$$\Re$$','$$\Im$$'); set(legend_handle,'Interpreter','latex');
subplot(3,3,5);
hold all
plot(f, real(em22),'r')
plot(f, imag(em22),'b')
hold off; xlabel(['Frequency [10^{' num2str(fpower) '} Hz]'])
grid on; ylabel('$$ \alpha_{\rm em}^{22} \rm{[m^2 s]}$$','Interpreter','latex')
legend_handle=legend('$$\Re$$','$$\Im$$'); set(legend_handle,'Interpreter','latex');
subplot(3,3,6);
hold all
plot(f, real(em23),'r')
plot(f, imag(em23),'b')
hold off; xlabel(['Frequency [10^{' num2str(fpower) '} Hz]'])
grid on; ylabel('$$ \alpha_{\rm em}^{23} \rm{[m^2 s]}$$','Interpreter','latex')
legend_handle=legend('$$\Re$$','$$\Im$$'); set(legend_handle,'Interpreter','latex');
subplot(3,3,7);
hold all
plot(f, real(em31),'r')
plot(f, imag(em31),'b')
hold off; xlabel(['Frequency [10^{' num2str(fpower) '} Hz]'])
grid on; ylabel('$$ \alpha_{\rm em}^{31} \rm{[m^2 s]}$$','Interpreter','latex')
legend_handle=legend('$$\Re$$','$$\Im$$'); set(legend_handle,'Interpreter','latex');
subplot(3,3,8);
hold all
plot(f, real(em32),'r')
plot(f, imag(em32),'b')
hold off; xlabel(['Frequency [10^{' num2str(fpower) '} Hz]'])
grid on; ylabel('$$ \alpha_{\rm em}^{32} \rm{[m^2 s]}$$','Interpreter','latex')
legend_handle=legend('$$\Re$$','$$\Im$$'); set(legend_handle,'Interpreter','latex');
subplot(3,3,9);
hold all
plot(f, real(em33),'r')
plot(f, imag(em33),'b')
hold off; xlabel(['Frequency [10^{' num2str(fpower) '} Hz]'])
grid on; ylabel('$$ \alpha_{\rm em}^{33} \rm{[m^2 s]}$$','Interpreter','latex')
legend_handle=legend('$$\Re$$','$$\Im$$'); set(legend_handle,'Interpreter','latex');

%% Plotting 9 components of alpha_me tensor

figure10 = figure('Name','9 components of alpha_me tensor');
set(0,'defaultlinelinewidth',3)
subplot(3,3,1);
hold all
plot(f, real(me11),'r')
plot(f, imag(me11),'b')
hold off; xlabel(['Frequency [10^{' num2str(fpower) '} Hz]'])
grid on; ylabel('$$ \alpha_{\rm me}^{11} \rm{[m^2 s]}$$','Interpreter','latex')
legend_handle=legend('$$\Re$$','$$\Im$$'); set(legend_handle,'Interpreter','latex');
subplot(3,3,2);
hold all
plot(f, real(me12),'r')
plot(f, imag(me12),'b')
hold off; xlabel(['Frequency [10^{' num2str(fpower) '} Hz]'])
grid on; ylabel('$$ \alpha_{\rm me}^{12} \rm{[m^2 s]}$$','Interpreter','latex')
legend_handle=legend('$$\Re$$','$$\Im$$'); set(legend_handle,'Interpreter','latex');
subplot(3,3,3);
hold all
plot(f, real(me13),'r')
plot(f, imag(me13),'b')
hold off; xlabel(['Frequency [10^{' num2str(fpower) '} Hz]'])
grid on; ylabel('$$ \alpha_{\rm me}^{13} \rm{[m^2 s]}$$','Interpreter','latex')
legend_handle=legend('$$\Re$$','$$\Im$$'); set(legend_handle,'Interpreter','latex');
subplot(3,3,4);
hold all
plot(f, real(me21),'r')
plot(f, imag(me21),'b')
hold off; xlabel(['Frequency [10^{' num2str(fpower) '} Hz]'])
grid on; ylabel('$$ \alpha_{\rm me}^{21} \rm{[m^2 s]}$$','Interpreter','latex')
legend_handle=legend('$$\Re$$','$$\Im$$'); set(legend_handle,'Interpreter','latex');
subplot(3,3,5);
hold all
plot(f, real(me22),'r')
plot(f, imag(me22),'b')
hold off; xlabel(['Frequency [10^{' num2str(fpower) '} Hz]'])
grid on; ylabel('$$ \alpha_{\rm me}^{22} \rm{[m^2 s]}$$','Interpreter','latex')
legend_handle=legend('$$\Re$$','$$\Im$$'); set(legend_handle,'Interpreter','latex');
subplot(3,3,6);
hold all
plot(f, real(me23),'r')
plot(f, imag(me23),'b')
hold off; xlabel(['Frequency [10^{' num2str(fpower) '} Hz]'])
grid on; ylabel('$$ \alpha_{\rm me}^{23} \rm{[m^2 s]}$$','Interpreter','latex')
legend_handle=legend('$$\Re$$','$$\Im$$'); set(legend_handle,'Interpreter','latex');
subplot(3,3,7);
hold all
plot(f, real(me31),'r')
plot(f, imag(me31),'b')
hold off; xlabel(['Frequency [10^{' num2str(fpower) '} Hz]'])
grid on; ylabel('$$ \alpha_{\rm me}^{31} \rm{[m^2 s]}$$','Interpreter','latex')
legend_handle=legend('$$\Re$$','$$\Im$$'); set(legend_handle,'Interpreter','latex');
subplot(3,3,8);
hold all
plot(f, real(me32),'r')
plot(f, imag(me32),'b')
hold off; xlabel(['Frequency [10^{' num2str(fpower) '} Hz]'])
grid on; ylabel('$$ \alpha_{\rm me}^{32} \rm{[m^2 s]}$$','Interpreter','latex')
legend_handle=legend('$$\Re$$','$$\Im$$'); set(legend_handle,'Interpreter','latex');
subplot(3,3,9);
hold all
plot(f, real(me33),'r')
plot(f, imag(me33),'b')
hold off; xlabel(['Frequency [10^{' num2str(fpower) '} Hz]'])
grid on; ylabel('$$ \alpha_{\rm me}^{33} \rm{[m^2 s]}$$','Interpreter','latex')
legend_handle=legend('$$\Re$$','$$\Im$$'); set(legend_handle,'Interpreter','latex');

%% Plotting 9 components of alpha_mm tensor

figure11 = figure('Name','9 components of alpha_mm tensor');
set(0,'defaultlinelinewidth',3)
subplot(3,3,1);
hold all
plot(f,1./eta.*real(mm11),'r')
plot(f,1./eta.*imag(mm11),'b')
hold off; xlabel(['Frequency [10^{' num2str(fpower) '} Hz]'])
grid on; ylabel('$$ 1/\eta_0 \alpha_{\rm mm}^{11} \rm{[m^2 s]}$$','Interpreter','latex')
legend_handle=legend('$$\Re$$','$$\Im$$'); set(legend_handle,'Interpreter','latex');
subplot(3,3,2);
hold all
plot(f,1./eta.*real(mm12),'r')
plot(f,1./eta.*imag(mm12),'b')
hold off; xlabel(['Frequency [10^{' num2str(fpower) '} Hz]'])
grid on; ylabel('$$ 1/\eta_0 \alpha_{\rm mm}^{12} \rm{[m^2 s]}$$','Interpreter','latex')
legend_handle=legend('$$\Re$$','$$\Im$$'); set(legend_handle,'Interpreter','latex');
subplot(3,3,3);
hold all
plot(f,1./eta.*real(mm13),'r')
plot(f,1./eta.*imag(mm13),'b')
hold off; xlabel(['Frequency [10^{' num2str(fpower) '} Hz]'])
grid on; ylabel('$$ 1/\eta_0 \alpha_{\rm mm}^{13} \rm{[m^2 s]}$$','Interpreter','latex')
legend_handle=legend('$$\Re$$','$$\Im$$'); set(legend_handle,'Interpreter','latex');
subplot(3,3,4);
hold all
plot(f,1./eta.*real(mm21),'r')
plot(f,1./eta.*imag(mm21),'b')
hold off; xlabel(['Frequency [10^{' num2str(fpower) '} Hz]'])
grid on; ylabel('$$ 1/\eta_0 \alpha_{\rm mm}^{21} \rm{[m^2 s]}$$','Interpreter','latex')
legend_handle=legend('$$\Re$$','$$\Im$$'); set(legend_handle,'Interpreter','latex');
subplot(3,3,5);
hold all
plot(f,1./eta.*real(mm22),'r')
plot(f,1./eta.*imag(mm22),'b')
hold off; xlabel(['Frequency [10^{' num2str(fpower) '} Hz]'])
grid on; ylabel('$$ 1/\eta_0 \alpha_{\rm mm}^{22} \rm{[m^2 s]}$$','Interpreter','latex')
legend_handle=legend('$$\Re$$','$$\Im$$'); set(legend_handle,'Interpreter','latex');
subplot(3,3,6);
hold all
plot(f,1./eta.*real(mm23),'r')
plot(f,1./eta.*imag(mm23),'b')
hold off; xlabel(['Frequency [10^{' num2str(fpower) '} Hz]'])
grid on; ylabel('$$ 1/\eta_0 \alpha_{\rm mm}^{23} \rm{[m^2 s]}$$','Interpreter','latex')
legend_handle=legend('$$\Re$$','$$\Im$$'); set(legend_handle,'Interpreter','latex');
subplot(3,3,7);
hold all
plot(f,1./eta.*real(mm31),'r')
plot(f,1./eta.*imag(mm31),'b')
hold off; xlabel(['Frequency [10^{' num2str(fpower) '} Hz]'])
grid on; ylabel('$$ 1/\eta_0 \alpha_{\rm mm}^{31} \rm{[m^2 s]}$$','Interpreter','latex')
legend_handle=legend('$$\Re$$','$$\Im$$'); set(legend_handle,'Interpreter','latex');
subplot(3,3,8);
hold all
plot(f,1./eta.*real(mm32),'r')
plot(f,1./eta.*imag(mm32),'b')
hold off; xlabel(['Frequency [10^{' num2str(fpower) '} Hz]'])
grid on; ylabel('$$ 1/\eta_0 \alpha_{\rm mm}^{32} \rm{[m^2 s]}$$','Interpreter','latex')
legend_handle=legend('$$\Re$$','$$\Im$$'); set(legend_handle,'Interpreter','latex');
subplot(3,3,9);
hold all
plot(f,1./eta.*real(mm33),'r')
plot(f,1./eta.*imag(mm33),'b')
hold off; xlabel(['Frequency [10^{' num2str(fpower) '} Hz]'])
grid on; ylabel('$$ 1/\eta_0 \alpha_{\rm mm}^{33} \rm{[m^2 s]}$$','Interpreter','latex')
legend_handle=legend('$$\Re$$','$$\Im$$'); set(legend_handle,'Interpreter','latex');

%% Decomposition to reciprocal and nonreciprocal parts
%NOTE: From now on the polarizabilities are normalized by impedance 
ee11r=eta.*0.5* 2*ee11;       ee22r=eta.*0.5* 2*ee22;       ee33r=eta.*0.5* 2*ee33; 
ee12r=eta.*0.5* (ee12+ee21);  ee13r=eta.*0.5* (ee13+ee31);  ee23r=eta.*0.5* (ee23+ee32); 
ee21r=ee12r;             ee31r=ee13r;             ee32r=ee23r; 

ee11n=eta.*zeros(1,length(ee11),'uint32');ee22n=eta.*zeros(1,length(ee11),'uint32'); ee33n=eta.*zeros(1,length(ee11),'uint32');
ee12n=eta.*0.5* (ee12-ee21);  ee13n=eta.*0.5* (ee13-ee31);  ee23n=eta.*0.5* (ee23-ee32); 
ee21n=-ee12n;                 ee31n=-ee13n;                 ee32n=-ee23n; 

mm11r=1./eta.*0.5* 2*mm11;       mm22r=1./eta.*0.5* 2*mm22;       mm33r=1./eta.*0.5* 2*mm33; 
mm12r=1./eta.*0.5* (mm12+mm21);  mm13r=1./eta.*0.5* (mm13+mm31);  mm23r=1./eta.*0.5* (mm23+mm32); 
mm21r=mm12r;             mm31r=mm13r;             mm32r=mm23r; 

mm11n=1./eta.*zeros(1,length(ee11),'uint32');  mm22n=1./eta.*zeros(1,length(ee11),'uint32');  mm33n=1./eta.*zeros(1,length(ee11),'uint32');
mm12n=1./eta.*0.5* (mm12-mm21);  mm13n=1./eta.*0.5* (mm13-mm31);  mm23n=1./eta.*0.5* (mm23-mm32); 
mm21n=-mm12n;            mm31n=-mm13n;            mm32n=-mm23n; 

em11r=0.5* (em11-me11);  em12r=0.5* (em12-me21);  em13r=0.5* (em13-me31);
em21r=0.5* (em21-me12);  em22r=0.5* (em22-me22);  em23r=0.5* (em23-me32); 
em31r=0.5* (em31-me13);  em32r=0.5* (em32-me23);  em33r=0.5* (em33-me33); 

em11n=0.5* (em11+me11);  em12n=0.5* (em12+me21);  em13n=0.5* (em13+me31);
em21n=0.5* (em21+me12);  em22n=0.5* (em22+me22);  em23n=0.5* (em23+me32); 
em31n=0.5* (em31+me13);  em32n=0.5* (em32+me23);  em33n=0.5* (em33+me33); 

% me11r=0.5* (me11-em11);  me12r=0.5* (me12-em21);  me13r=0.5* (me13-em31);
% me21r=0.5* (me21-em12);  me22r=0.5* (me22-em22);  me23r=0.5* (me23-em32); 
% me31r=0.5* (me31-em13);  me32r=0.5* (me32-em23);  me33r=0.5* (me33-em33); 
% 
% me11n=0.5* (me11+em11);  me12n=0.5* (me12+em21);  me13n=0.5* (me13+em31);
% me21n=0.5* (me21+em12);  me22n=0.5* (me22+em22);  me23n=0.5* (me23+em32); 
% me31n=0.5* (me31+em13);  me32n=0.5* (me32+em23);  me33n=0.5* (me33+em33); 

%% Plotting alpha_ee, reciprocal and nonreciprocal parts

figure12 = figure('Name','Alpha_ee, reciprocal and nonreciprocal parts');
set(0,'defaultlinelinewidth',3)
subplot(3,3,1);
hold all
plot(f, real(ee11r),'r')
plot(f, imag(ee11r),'b')
plot(f, real(ee11n),'r:')
plot(f, imag(ee11n),'b:')
hold off; xlabel(['Frequency [10^{' num2str(fpower) '} Hz]'])
grid on; ylabel('$$ \eta_0 \alpha_{\rm ee}^{11} \rm{[m^2 s]}$$','Interpreter','latex')
legend_handle=legend('$$\Re {\rm (rec.)}$$','$$\Im {\rm (rec.)}$$','$$\Re {\rm (nonr.)}$$','$$\Im {\rm (nonr.)}$$'); 
set(legend_handle,'Interpreter','latex','location','northoutside','orientation','horizontal','NumColumns',2,'FontSize',9);
subplot(3,3,2);
hold all
plot(f, real(ee12r),'r')
plot(f, imag(ee12r),'b')
plot(f, real(ee12n),'r:')
plot(f, imag(ee12n),'b:')
hold off; xlabel(['Frequency [10^{' num2str(fpower) '} Hz]'])
grid on; ylabel('$$ \eta_0 \alpha_{\rm ee}^{12} \rm{[m^2 s]}$$','Interpreter','latex')
legend_handle=legend('$$\Re {\rm (rec.)}$$','$$\Im {\rm (rec.)}$$','$$\Re {\rm (nonr.)}$$','$$\Im {\rm (nonr.)}$$'); 
set(legend_handle,'Interpreter','latex','location','northoutside','orientation','horizontal','NumColumns',2,'FontSize',9);

subplot(3,3,3);
hold all
plot(f, real(ee13r),'r')
plot(f, imag(ee13r),'b')
plot(f, real(ee13n),'r:')
plot(f, imag(ee13n),'b:')
hold off; xlabel(['Frequency [10^{' num2str(fpower) '} Hz]'])
grid on; ylabel('$$ \eta_0 \alpha_{\rm ee}^{13} \rm{[m^2 s]}$$','Interpreter','latex')
legend_handle=legend('$$\Re {\rm (rec.)}$$','$$\Im {\rm (rec.)}$$','$$\Re {\rm (nonr.)}$$','$$\Im {\rm (nonr.)}$$'); 
set(legend_handle,'Interpreter','latex','location','northoutside','orientation','horizontal','NumColumns',2,'FontSize',9);

subplot(3,3,4);
hold all
plot(f, real(ee21r),'r')
plot(f, imag(ee21r),'b')
plot(f, real(ee21n),'r:')
plot(f, imag(ee21n),'b:')
hold off; xlabel(['Frequency [10^{' num2str(fpower) '} Hz]'])
grid on; ylabel('$$ \eta_0 \alpha_{\rm ee}^{21} \rm{[m^2 s]}$$','Interpreter','latex')
legend_handle=legend('$$\Re {\rm (rec.)}$$','$$\Im {\rm (rec.)}$$','$$\Re {\rm (nonr.)}$$','$$\Im {\rm (nonr.)}$$'); 
set(legend_handle,'Interpreter','latex','location','northoutside','orientation','horizontal','NumColumns',2,'FontSize',9);

subplot(3,3,5);
hold all
plot(f, real(ee22r),'r')
plot(f, imag(ee22r),'b')
plot(f, real(ee22n),'r:')
plot(f, imag(ee22n),'b:')
hold off; xlabel(['Frequency [10^{' num2str(fpower) '} Hz]'])
grid on; ylabel('$$ \eta_0 \alpha_{\rm ee}^{22} \rm{[m^2 s]}$$','Interpreter','latex')
legend_handle=legend('$$\Re {\rm (rec.)}$$','$$\Im {\rm (rec.)}$$','$$\Re {\rm (nonr.)}$$','$$\Im {\rm (nonr.)}$$'); 
set(legend_handle,'Interpreter','latex','location','northoutside','orientation','horizontal','NumColumns',2,'FontSize',9);

subplot(3,3,6);
hold all
plot(f, real(ee23r),'r')
plot(f, imag(ee23r),'b')
plot(f, real(ee23n),'r:')
plot(f, imag(ee23n),'b:')
hold off; xlabel(['Frequency [10^{' num2str(fpower) '} Hz]'])
grid on; ylabel('$$ \eta_0 \alpha_{\rm ee}^{23} \rm{[m^2 s]}$$','Interpreter','latex')
legend_handle=legend('$$\Re {\rm (rec.)}$$','$$\Im {\rm (rec.)}$$','$$\Re {\rm (nonr.)}$$','$$\Im {\rm (nonr.)}$$'); 
set(legend_handle,'Interpreter','latex','location','northoutside','orientation','horizontal','NumColumns',2,'FontSize',9);

subplot(3,3,7);
hold all
plot(f, real(ee31r),'r')
plot(f, imag(ee31r),'b')
plot(f, real(ee31n),'r:')
plot(f, imag(ee31n),'b:')
hold off; xlabel(['Frequency [10^{' num2str(fpower) '} Hz]'])
grid on; ylabel('$$ \eta_0 \alpha_{\rm ee}^{31} \rm{[m^2 s]}$$','Interpreter','latex')
legend_handle=legend('$$\Re {\rm (rec.)}$$','$$\Im {\rm (rec.)}$$','$$\Re {\rm (nonr.)}$$','$$\Im {\rm (nonr.)}$$'); 
set(legend_handle,'Interpreter','latex','location','northoutside','orientation','horizontal','NumColumns',2,'FontSize',9);

subplot(3,3,8);
hold all
plot(f, real(ee32r),'r')
plot(f, imag(ee32r),'b')
plot(f, real(ee32n),'r:')
plot(f, imag(ee32n),'b:')
hold off; xlabel(['Frequency [10^{' num2str(fpower) '} Hz]'])
grid on; ylabel('$$ \eta_0 \alpha_{\rm ee}^{32} \rm{[m^2 s]}$$','Interpreter','latex')
legend_handle=legend('$$\Re {\rm (rec.)}$$','$$\Im {\rm (rec.)}$$','$$\Re {\rm (nonr.)}$$','$$\Im {\rm (nonr.)}$$'); 
set(legend_handle,'Interpreter','latex','location','northoutside','orientation','horizontal','NumColumns',2,'FontSize',9);

subplot(3,3,9);
hold all
plot(f, real(ee33r),'r')
plot(f, imag(ee33r),'b')
plot(f, real(ee33n),'r:')
plot(f, imag(ee33n),'b:')
hold off; xlabel(['Frequency [10^{' num2str(fpower) '} Hz]'])
grid on; ylabel('$$ \eta_0 \alpha_{\rm ee}^{33} \rm{[m^2 s]}$$','Interpreter','latex')
legend_handle=legend('$$\Re {\rm (rec.)}$$','$$\Im {\rm (rec.)}$$','$$\Re {\rm (nonr.)}$$','$$\Im {\rm (nonr.)}$$'); 
set(legend_handle,'Interpreter','latex','location','northoutside','orientation','horizontal','NumColumns',2,'FontSize',9);

%% Plotting alpha_mm, reciprocal and nonreciprocal parts

figure13 = figure('Name','Alpha_mm, reciprocal and nonreciprocal parts');
set(0,'defaultlinelinewidth',3)
subplot(3,3,1);
hold all
plot(f, real(mm11r),'r')
plot(f, imag(mm11r),'b')
plot(f, real(mm11n),'r:')
plot(f, imag(mm11n),'b:')
hold off; xlabel(['Frequency [10^{' num2str(fpower) '} Hz]'])
grid on; ylabel('$$ 1/\eta_0 \alpha_{\rm mm}^{11} \rm{[m^2 s]}$$','Interpreter','latex')
legend_handle=legend('$$\Re {\rm (rec.)}$$','$$\Im {\rm (rec.)}$$','$$\Re {\rm (nonr.)}$$','$$\Im {\rm (nonr.)}$$'); 
set(legend_handle,'Interpreter','latex','location','northoutside','orientation','horizontal','NumColumns',2,'FontSize',9);

subplot(3,3,2);
hold all
plot(f, real(mm12r),'r')
plot(f, imag(mm12r),'b')
plot(f, real(mm12n),'r:')
plot(f, imag(mm12n),'b:')
hold off; xlabel(['Frequency [10^{' num2str(fpower) '} Hz]'])
grid on; ylabel('$$ 1/\eta_0 \alpha_{\rm mm}^{12} \rm{[m^2 s]}$$','Interpreter','latex')
legend_handle=legend('$$\Re {\rm (rec.)}$$','$$\Im {\rm (rec.)}$$','$$\Re {\rm (nonr.)}$$','$$\Im {\rm (nonr.)}$$'); 
set(legend_handle,'Interpreter','latex','location','northoutside','orientation','horizontal','NumColumns',2,'FontSize',9);

subplot(3,3,3);
hold all
plot(f, real(mm13r),'r')
plot(f, imag(mm13r),'b')
plot(f, real(mm13n),'r:')
plot(f, imag(mm13n),'b:')
hold off; xlabel(['Frequency [10^{' num2str(fpower) '} Hz]'])
grid on; ylabel('$$ 1/\eta_0 \alpha_{\rm mm}^{13} \rm{[m^2 s]}$$','Interpreter','latex')
legend_handle=legend('$$\Re {\rm (rec.)}$$','$$\Im {\rm (rec.)}$$','$$\Re {\rm (nonr.)}$$','$$\Im {\rm (nonr.)}$$'); 
set(legend_handle,'Interpreter','latex','location','northoutside','orientation','horizontal','NumColumns',2,'FontSize',9);

subplot(3,3,4);
hold all
plot(f, real(mm21r),'r')
plot(f, imag(mm21r),'b')
plot(f, real(mm21n),'r:')
plot(f, imag(mm21n),'b:')
hold off; xlabel(['Frequency [10^{' num2str(fpower) '} Hz]'])
grid on; ylabel('$$ 1/\eta_0 \alpha_{\rm mm}^{21} \rm{[m^2 s]}$$','Interpreter','latex')
legend_handle=legend('$$\Re {\rm (rec.)}$$','$$\Im {\rm (rec.)}$$','$$\Re {\rm (nonr.)}$$','$$\Im {\rm (nonr.)}$$'); 
set(legend_handle,'Interpreter','latex','location','northoutside','orientation','horizontal','NumColumns',2,'FontSize',9);

subplot(3,3,5);
hold all
plot(f, real(mm22r),'r')
plot(f, imag(mm22r),'b')
plot(f, real(mm22n),'r:')
plot(f, imag(mm22n),'b:')
hold off; xlabel(['Frequency [10^{' num2str(fpower) '} Hz]'])
grid on; ylabel('$$ 1/\eta_0 \alpha_{\rm mm}^{22} \rm{[m^2 s]}$$','Interpreter','latex')
legend_handle=legend('$$\Re {\rm (rec.)}$$','$$\Im {\rm (rec.)}$$','$$\Re {\rm (nonr.)}$$','$$\Im {\rm (nonr.)}$$'); 
set(legend_handle,'Interpreter','latex','location','northoutside','orientation','horizontal','NumColumns',2,'FontSize',9);

subplot(3,3,6);
hold all
plot(f, real(mm23r),'r')
plot(f, imag(mm23r),'b')
plot(f, real(mm23n),'r:')
plot(f, imag(mm23n),'b:')
hold off; xlabel(['Frequency [10^{' num2str(fpower) '} Hz]'])
grid on; ylabel('$$ 1/\eta_0 \alpha_{\rm mm}^{23} \rm{[m^2 s]}$$','Interpreter','latex')
legend_handle=legend('$$\Re {\rm (rec.)}$$','$$\Im {\rm (rec.)}$$','$$\Re {\rm (nonr.)}$$','$$\Im {\rm (nonr.)}$$'); 
set(legend_handle,'Interpreter','latex','location','northoutside','orientation','horizontal','NumColumns',2,'FontSize',9);

subplot(3,3,7);
hold all
plot(f, real(mm31r),'r')
plot(f, imag(mm31r),'b')
plot(f, real(mm31n),'r:')
plot(f, imag(mm31n),'b:')
hold off; xlabel(['Frequency [10^{' num2str(fpower) '} Hz]'])
grid on; ylabel('$$ 1/\eta_0 \alpha_{\rm mm}^{31} \rm{[m^2 s]}$$','Interpreter','latex')
legend_handle=legend('$$\Re {\rm (rec.)}$$','$$\Im {\rm (rec.)}$$','$$\Re {\rm (nonr.)}$$','$$\Im {\rm (nonr.)}$$'); 
set(legend_handle,'Interpreter','latex','location','northoutside','orientation','horizontal','NumColumns',2,'FontSize',9);

subplot(3,3,8);
hold all
plot(f, real(mm32r),'r')
plot(f, imag(mm32r),'b')
plot(f, real(mm32n),'r:')
plot(f, imag(mm32n),'b:')
hold off; xlabel(['Frequency [10^{' num2str(fpower) '} Hz]'])
grid on; ylabel('$$ 1/\eta_0 \alpha_{\rm mm}^{32} \rm{[m^2 s]}$$','Interpreter','latex')
legend_handle=legend('$$\Re {\rm (rec.)}$$','$$\Im {\rm (rec.)}$$','$$\Re {\rm (nonr.)}$$','$$\Im {\rm (nonr.)}$$'); 
set(legend_handle,'Interpreter','latex','location','northoutside','orientation','horizontal','NumColumns',2,'FontSize',9);

subplot(3,3,9);
hold all
plot(f, real(mm33r),'r')
plot(f, imag(mm33r),'b')
plot(f, real(mm33n),'r:')
plot(f, imag(mm33n),'b:')
hold off; xlabel(['Frequency [10^{' num2str(fpower) '} Hz]'])
grid on; ylabel('$$ 1/\eta_0 \alpha_{\rm mm}^{33} \rm{[m^2 s]}$$','Interpreter','latex')
legend_handle=legend('$$\Re {\rm (rec.)}$$','$$\Im {\rm (rec.)}$$','$$\Re {\rm (nonr.)}$$','$$\Im {\rm (nonr.)}$$'); 
set(legend_handle,'Interpreter','latex','location','northoutside','orientation','horizontal','NumColumns',2,'FontSize',9);

%% Plotting alpha_em, reciprocal and nonreciprocal parts
figure14 = figure('Name','Alpha_em, reciprocal and nonreciprocal parts');
set(0,'defaultlinelinewidth',3)
subplot(3,3,1);
hold all
plot(f, real(em11r),'r')
plot(f, imag(em11r),'b')
plot(f, real(em11n),'r:')
plot(f, imag(em11n),'b:')
hold off; xlabel(['Frequency [10^{' num2str(fpower) '} Hz]'])
grid on; ylabel('$$ \alpha_{\rm em}^{11} \rm{[m^2 s]}$$','Interpreter','latex')
legend_handle=legend('$$\Re {\rm (rec.)}$$','$$\Im {\rm (rec.)}$$','$$\Re {\rm (nonr.)}$$','$$\Im {\rm (nonr.)}$$'); 
set(legend_handle,'Interpreter','latex','location','northoutside','orientation','horizontal','NumColumns',2,'FontSize',9);

subplot(3,3,2);
hold all
plot(f, real(em12r),'r')
plot(f, imag(em12r),'b')
plot(f, real(em12n),'r:')
plot(f, imag(em12n),'b:')
hold off; xlabel(['Frequency [10^{' num2str(fpower) '} Hz]'])
grid on; ylabel('$$ \alpha_{\rm em}^{12} \rm{[m^2 s]}$$','Interpreter','latex')
legend_handle=legend('$$\Re {\rm (rec.)}$$','$$\Im {\rm (rec.)}$$','$$\Re {\rm (nonr.)}$$','$$\Im {\rm (nonr.)}$$'); 
set(legend_handle,'Interpreter','latex','location','northoutside','orientation','horizontal','NumColumns',2,'FontSize',9);

subplot(3,3,3);
hold all
plot(f, real(em13r),'r')
plot(f, imag(em13r),'b')
plot(f, real(em13n),'r:')
plot(f, imag(em13n),'b:')
hold off; xlabel(['Frequency [10^{' num2str(fpower) '} Hz]'])
grid on; ylabel('$$ \alpha_{\rm em}^{13} \rm{[m^2 s]}$$','Interpreter','latex')
legend_handle=legend('$$\Re {\rm (rec.)}$$','$$\Im {\rm (rec.)}$$','$$\Re {\rm (nonr.)}$$','$$\Im {\rm (nonr.)}$$'); 
set(legend_handle,'Interpreter','latex','location','northoutside','orientation','horizontal','NumColumns',2,'FontSize',9);

subplot(3,3,4);
hold all
plot(f, real(em21r),'r')
plot(f, imag(em21r),'b')
plot(f, real(em21n),'r:')
plot(f, imag(em21n),'b:')
hold off; xlabel(['Frequency [10^{' num2str(fpower) '} Hz]'])
grid on; ylabel('$$ \alpha_{\rm em}^{21} \rm{[m^2 s]}$$','Interpreter','latex')
legend_handle=legend('$$\Re {\rm (rec.)}$$','$$\Im {\rm (rec.)}$$','$$\Re {\rm (nonr.)}$$','$$\Im {\rm (nonr.)}$$'); 
set(legend_handle,'Interpreter','latex','location','northoutside','orientation','horizontal','NumColumns',2,'FontSize',9);

subplot(3,3,5);
hold all
plot(f, real(em22r),'r')
plot(f, imag(em22r),'b')
plot(f, real(em22n),'r:')
plot(f, imag(em22n),'b:')
hold off; xlabel(['Frequency [10^{' num2str(fpower) '} Hz]'])
grid on; ylabel('$$ \alpha_{\rm em}^{22} \rm{[m^2 s]}$$','Interpreter','latex')
legend_handle=legend('$$\Re {\rm (rec.)}$$','$$\Im {\rm (rec.)}$$','$$\Re {\rm (nonr.)}$$','$$\Im {\rm (nonr.)}$$'); 
set(legend_handle,'Interpreter','latex','location','northoutside','orientation','horizontal','NumColumns',2,'FontSize',9);

subplot(3,3,6);
hold all
plot(f, real(em23r),'r')
plot(f, imag(em23r),'b')
plot(f, real(em23n),'r:')
plot(f, imag(em23n),'b:')
hold off; xlabel(['Frequency [10^{' num2str(fpower) '} Hz]'])
grid on; ylabel('$$ \alpha_{\rm em}^{23} \rm{[m^2 s]}$$','Interpreter','latex')
legend_handle=legend('$$\Re {\rm (rec.)}$$','$$\Im {\rm (rec.)}$$','$$\Re {\rm (nonr.)}$$','$$\Im {\rm (nonr.)}$$'); 
set(legend_handle,'Interpreter','latex','location','northoutside','orientation','horizontal','NumColumns',2,'FontSize',9);

subplot(3,3,7);
hold all
plot(f, real(em31r),'r')
plot(f, imag(em31r),'b')
plot(f, real(em31n),'r:')
plot(f, imag(em31n),'b:')
hold off; xlabel(['Frequency [10^{' num2str(fpower) '} Hz]'])
grid on; ylabel('$$ \alpha_{\rm em}^{31} \rm{[m^2 s]}$$','Interpreter','latex')
legend_handle=legend('$$\Re {\rm (rec.)}$$','$$\Im {\rm (rec.)}$$','$$\Re {\rm (nonr.)}$$','$$\Im {\rm (nonr.)}$$'); 
set(legend_handle,'Interpreter','latex','location','northoutside','orientation','horizontal','NumColumns',2,'FontSize',9);

subplot(3,3,8);
hold all
plot(f, real(em32r),'r')
plot(f, imag(em32r),'b')
plot(f, real(em32n),'r:')
plot(f, imag(em32n),'b:')
hold off; xlabel(['Frequency [10^{' num2str(fpower) '} Hz]'])
grid on; ylabel('$$ \alpha_{\rm em}^{32} \rm{[m^2 s]}$$','Interpreter','latex')
legend_handle=legend('$$\Re {\rm (rec.)}$$','$$\Im {\rm (rec.)}$$','$$\Re {\rm (nonr.)}$$','$$\Im {\rm (nonr.)}$$'); 
set(legend_handle,'Interpreter','latex','location','northoutside','orientation','horizontal','NumColumns',2,'FontSize',9);

subplot(3,3,9);
hold all
plot(f, real(em33r),'r')
plot(f, imag(em33r),'b')
plot(f, real(em33n),'r:')
plot(f, imag(em33n),'b:')
hold off; xlabel(['Frequency [10^{' num2str(fpower) '} Hz]'])
grid on; ylabel('$$ \alpha_{\rm em}^{33} \rm{[m^2 s]}$$','Interpreter','latex')
legend_handle=legend('$$\Re {\rm (rec.)}$$','$$\Im {\rm (rec.)}$$','$$\Re {\rm (nonr.)}$$','$$\Im {\rm (nonr.)}$$'); 
set(legend_handle,'Interpreter','latex','location','northoutside','orientation','horizontal','NumColumns',2,'FontSize',9);

%% Plotting strongest polarizability components in descending order
% 1     2     3     4     5     6     7     8     9         10    11    12    13    14    15    16    17    18 
% ee11r ee22r ee33r ee12r ee13r ee23r ee12n ee13n ee23n     mm11r mm22r mm33r mm12r mm13r mm23r mm12n mm13n mm23n
% 19    20    21    22    23    24    25    26    27        28    29    30    31    32    33    34    35    36   
% em11r em12r em13r em21r em22r em23r em31r em32r em33r     em11n em12n em13n em21n em22n em23n em31n em32n em33n
names = ["$$\eta_0\alpha_{\rm ee, rec}^{11}$$","$$\eta_0\alpha_{\rm ee, rec}^{22}$$","$$\eta_0\alpha_{\rm ee, rec}^{33}$$",...
    "$$\eta_0\alpha_{\rm ee, rec}^{12}$$","$$\eta_0\alpha_{\rm ee, rec}^{13}$$","$$\eta_0\alpha_{\rm ee, rec}^{23}$$",...
    "$$\eta_0\alpha_{\rm ee, nonr}^{12}$$","$$\eta_0\alpha_{\rm ee, nonr}^{13}$$","$$\eta_0\alpha_{\rm ee, nonr}^{23}$$",...
    "$$1/\eta_0\alpha_{\rm mm, rec}^{11}$$","$$1/\eta_0\alpha_{\rm mm, rec}^{22}$$","$$1/\eta_0\alpha_{\rm mm, rec}^{33}$$",...
    "$$1/\eta_0\alpha_{\rm mm, rec}^{12}$$","$$1/\eta_0\alpha_{\rm mm, rec}^{13}$$","$$1/\eta_0\alpha_{\rm mm, rec}^{23}$$",...
    "$$1/\eta_0\alpha_{\rm mm, nonr}^{12}$$","$$1/\eta_0\alpha_{\rm mm, nonr}^{13}$$","$$1/\eta_0\alpha_{\rm mm, nonr}^{23}$$",...
    "$$\alpha_{\rm em, rec}^{11}$$","$$\alpha_{\rm em, rec}^{12}$$","$$\alpha_{\rm em, rec}^{13}$$",...
    "$$\alpha_{\rm em, rec}^{21}$$","$$\alpha_{\rm em, rec}^{22}$$","$$\alpha_{\rm em, rec}^{23}$$",...
    "$$\alpha_{\rm em, rec}^{31}$$","$$\alpha_{\rm em, rec}^{32}$$","$$\alpha_{\rm em, rec}^{33}$$",...
    "$$\alpha_{\rm em, nonr}^{11}$$","$$\alpha_{\rm em, nonr}^{12}$$","$$\alpha_{\rm em, nonr}^{13}$$",...
    "$$\alpha_{\rm em, nonr}^{21}$$","$$\alpha_{\rm em, nonr}^{22}$$","$$\alpha_{\rm em, nonr}^{23}$$",...
    "$$\alpha_{\rm em, nonr}^{31}$$","$$\alpha_{\rm em, nonr}^{32}$$","$$\alpha_{\rm em, nonr}^{33}$$"];

Alpha=[ee11r,ee22r,ee33r,ee12r,ee13r,ee23r,ee12n,ee13n,ee23n,mm11r,mm22r,mm33r,mm12r,mm13r,mm23r,mm12n,mm13n,mm23n,...
  em11r,em12r,em13r,em21r,em22r,em23r,em31r,em32r,em33r,em11n,em12n,em13n,em21n,em22n,em23n,em31n,em32n,em33n];
maxAlpha=max(abs(Alpha));
[B,I] = sort(maxAlpha,'descend');

figure15 = figure('Name','First 9 components sorted by magnitude in descending order');
set(0,'defaultlinelinewidth',3)
set(0,'DefaultAxesFontSize', 10)
for ii=1:9
subplot(3,3,ii);
plot(f,real(Alpha(:,I(ii))),'r')
hold on
plot(f,imag(Alpha(:,I(ii))),'b')
grid on;
ylabel(names(I(ii)),'FontSize',14,'Interpreter','latex')
title(['$$|$$'+names(I(ii))+'$$|_{\rm max} = $$'+num2str(maxAlpha(I(ii)))],'FontSize',11,'Interpreter','latex');
xlabel(['Frequency [10^{' num2str(fpower) '} Hz]'])
legend_handle=legend('$$\Re$$','$$\Im$$'); 
set(legend_handle,'Interpreter','latex');
end;
 
%% Determining the relative error of the polarizability extraction

max_ampl=max(abs(Alpha(:,I(1))));

error_ee11=eta.*max(abs(ee11-ee11copy))/max_ampl;
error_ee22=eta.*max(abs(ee22-ee22copy))/max_ampl;
error_ee33=eta.*max(abs(ee33-ee33copy))/max_ampl;
error_em11=max(abs(em11-em11copy))/max_ampl;
error_em22=max(abs(em22-em22copy))/max_ampl;
error_em33=max(abs(em33-em33copy))/max_ampl;
error_me11=max(abs(me11-me11copy))/max_ampl;
error_me22=max(abs(me22-me22copy))/max_ampl;
error_me33=max(abs(me33-me33copy))/max_ampl;
error_mm11=1./eta.*max(abs(mm11-mm11copy))/max_ampl;
error_mm22=1./eta.*max(abs(mm22-mm22copy))/max_ampl;
error_mm33=1./eta.*max(abs(mm33-mm33copy))/max_ampl;

Alpha_f0=[ee11r(f0),ee22r(f0),ee33r(f0),ee12r(f0),ee13r(f0),ee23r(f0),ee12n(f0),ee13n(f0),ee23n(f0),mm11r(f0),mm22r(f0),mm33r(f0),mm12r(f0),mm13r(f0),mm23r(f0),mm12n(f0),mm13n(f0),mm23n(f0),...
  em11r(f0),em12r(f0),em13r(f0),em21r(f0),em22r(f0),em23r(f0),em31r(f0),em32r(f0),em33r(f0),em11n(f0),em12n(f0),em13n(f0),em21n(f0),em22n(f0),em23n(f0),em31n(f0),em32n(f0),em33n(f0)];
max_ampl_f0=max(max(abs(Alpha)));

error_ee11_f0=eta.*(abs(ee11(f0)-ee11copy(f0)))/max_ampl_f0;
error_ee22_f0=eta.*(abs(ee22(f0)-ee22copy(f0)))/max_ampl_f0;
error_ee33_f0=eta.*(abs(ee33(f0)-ee33copy(f0)))/max_ampl_f0;
error_em11_f0=(abs(em11(f0)-em11copy(f0)))/max_ampl_f0;
error_em22_f0=(abs(em22(f0)-em22copy(f0)))/max_ampl_f0;
error_em33_f0=(abs(em33(f0)-em33copy(f0)))/max_ampl_f0;
error_me11_f0=(abs(me11(f0)-me11copy(f0)))/max_ampl_f0;
error_me22_f0=(abs(me22(f0)-me22copy(f0)))/max_ampl_f0;
error_me33_f0=(abs(me33(f0)-me33copy(f0)))/max_ampl_f0;
error_mm11_f0=1./eta.*(abs(mm11(f0)-mm11copy(f0)))/max_ampl_f0;
error_mm22_f0=1./eta.*(abs(mm22(f0)-mm22copy(f0)))/max_ampl_f0;
error_mm33_f0=1./eta.*(abs(mm33(f0)-mm33copy(f0)))/max_ampl_f0;

'Relative error of the polarizability extraction for your meta-atom in the entire frequency range:'
relative_error=max([error_ee11 error_ee22 error_ee33 error_em11 error_em22 error_em33 ...
error_me11 error_me22 error_me33 error_mm11 error_mm22 error_mm33])

'Relative error of the polarizability extraction for your meta-atom at the chosen frequency:'
relative_error_f0=max([error_ee11_f0 error_ee22_f0 error_ee33_f0 error_em11_f0 error_em22_f0 error_em33_f0 ...
error_me11_f0 error_me22_f0 error_me33_f0 error_mm11_f0 error_mm22_f0 error_mm33_f0])

%% Modulus decomposition for alpha_ee reciprocal
% the tensor is symmetric and therefore is decomposed to diagonal and symmetric parts as: alpha=DeerI + S_eer 

% diagonal part is realized as 3 equal electric dipoles oriented along the
% x, y, and z axes. Their amplitudes are equal to Deer
Deer=(ee11r(f0)+ee22r(f0)+ee33r(f0))/3;

% symmetric part is realized as 6 electric dipoles (with different complex
% aimplitudes Seer_i) oriented along 6 real vectors veer_i
sr(1,1)=(2*ee11r(f0)-ee22r(f0)-ee33r(f0))/3;
sr(2,2)=(2*ee22r(f0)-ee11r(f0)-ee33r(f0))/3;
sr(3,3)=(2*ee33r(f0)-ee11r(f0)-ee22r(f0))/3;
sr(1,2)=ee12r(f0); sr(1,3)=ee13r(f0); sr(2,3)=ee23r(f0);
veer(:,1)=[1 0 0]'; veer(:,2)=[0 1 0]'; veer(:,3)=[0 0 1]';
for ii=1:2       %this cycle provides the best decomposition (with minimum number of inclusions)
    for jj=1:2
        for kk=1:2
            seer(ii,jj,kk)=abs(sr(1,1)-(-1)^(ii+1)*sr(1,2)-(-1)^(jj+1)*sr(1,3))+abs(sr(2,2)-(-1)^(ii+1)*sr(1,2)...
                -(-1)^(kk+1)*sr(2,3))+abs(sr(3,3)-(-1)^(jj+1)*sr(1,3)-(-1)^(kk+1)*sr(2,3));
        end;
    end
end;
[Meer,Ieer] = min(seer(:));
[iimin, jjmin, kkmin] = ind2sub(size(seer),Ieer);
Seer(1)=sr(1,1)-(-1)^(iimin+1)*sr(1,2)-(-1)^(jjmin+1)*sr(1,3);
Seer(2)=sr(2,2)-(-1)^(iimin+1)*sr(1,2)-(-1)^(kkmin+1)*sr(2,3);
Seer(3)=sr(3,3)-(-1)^(jjmin+1)*sr(1,3)-(-1)^(kkmin+1)*sr(2,3);
Seer(4)=2*(-1)^(iimin+1)*sr(1,2);
Seer(5)=2*(-1)^(jjmin+1)*sr(1,3);
Seer(6)=2*(-1)^(kkmin+1)*sr(2,3);
veer(:,4)=1/sqrt(2)*[1 (-1)^(iimin+1) 0]'; veer(:,5)=1/sqrt(2)*[1 0 (-1)^(jjmin+1)]'; veer(:,6)=1/sqrt(2)*[0 1 (-1)^(kkmin+1)]';

%% Modulus decomposition for alpha_ee non-reciprocal
% the tensor is antisymmetric and therefore is represented by only antisymmetric part alpha=A_een 

% the antisymmetric part is realized as 2 precessing electric dipoles: 
% one with purely real Aeen(1) and another with purely imaginary Aeen(2) 
% polarizability. The precessing axes are given by real vectors ween_i.

an(1,2)=ee12n(f0); an(1,3)=ee13n(f0); an(2,3)=ee23n(f0);
Aeen(1)=sqrt( (real(an(2,3)) )^2 + (real(an(1,3)) )^2 + (real(an(1,2)) )^2 );
Aeen(2)=1i*sqrt( (imag(an(2,3)) )^2 + (imag(an(1,3)) )^2 + (imag(an(1,2)) )^2 );
ween(:,1)=1/Aeen(1)*[-real(an(2,3))  real(an(1,3))  -real(an(1,2))]';
ween(:,2)=1/imag(Aeen(2))*[-imag(an(2,3))  imag(an(1,3))  -imag(an(1,2))]';

%% Modulus decomposition for alpha_mm reciprocal

% the tensor is symmetric and therefore is decomposed to diagonal and symmetric parts as: alpha=DmmrI + S_mmr 

% diagonal part is realized as 3 equal magnetic dipoles (or DSSR) oriented along the
% x, y, and z axes. Their amplitudes are equal to Dmmr
Dmmr=(mm11r(f0)+mm22r(f0)+mm33r(f0))/3;

% symmetric part is realized as 6 magnetic dipoles or DSRR (with different complex
% aimplitudes Smmr_i) oriented along 6 real vectors vmmr_i
sr(1,1)=(2*mm11r(f0)-mm22r(f0)-mm33r(f0))/3;
sr(2,2)=(2*mm22r(f0)-mm11r(f0)-mm33r(f0))/3;
sr(3,3)=(2*mm33r(f0)-mm11r(f0)-mm22r(f0))/3;
sr(1,2)=mm12r(f0); sr(1,3)=mm13r(f0); sr(2,3)=mm23r(f0);
vmmr(:,1)=[1 0 0]'; vmmr(:,2)=[0 1 0]'; vmmr(:,3)=[0 0 1]';
for ii=1:2       %this cycle provides the best decomposition (with minimum number of inclusions)
    for jj=1:2
        for kk=1:2
            smmr(ii,jj,kk)=abs(sr(1,1)-(-1)^(ii+1)*sr(1,2)-(-1)^(jj+1)*sr(1,3))+abs(sr(2,2)-(-1)^(ii+1)*sr(1,2)...
                -(-1)^(kk+1)*sr(2,3))+abs(sr(3,3)-(-1)^(jj+1)*sr(1,3)-(-1)^(kk+1)*sr(2,3));
        end;
    end
end;
[Mmmr,Immr] = min(smmr(:));
[iimin, jjmin, kkmin] = ind2sub(size(smmr),Immr);
Smmr(1)=sr(1,1)-(-1)^(iimin+1)*sr(1,2)-(-1)^(jjmin+1)*sr(1,3);
Smmr(2)=sr(2,2)-(-1)^(iimin+1)*sr(1,2)-(-1)^(kkmin+1)*sr(2,3);
Smmr(3)=sr(3,3)-(-1)^(jjmin+1)*sr(1,3)-(-1)^(kkmin+1)*sr(2,3);
Smmr(4)=2*(-1)^(iimin+1)*sr(1,2);
Smmr(5)=2*(-1)^(jjmin+1)*sr(1,3);
Smmr(6)=2*(-1)^(kkmin+1)*sr(2,3);
vmmr(:,4)=1/sqrt(2)*[1 (-1)^(iimin+1) 0]'; vmmr(:,5)=1/sqrt(2)*[1 0 (-1)^(jjmin+1)]'; vmmr(:,6)=1/sqrt(2)*[0 1 (-1)^(kkmin+1)]';

%% Modulus decomposition for alpha_mm non-reciprocal
% the tensor is antisymmetric and therefore is represented by only antisymmetric part alpha=A_mmn 

% the antisymmetric part is realized as 2 precessing magnetic dipoles: 
% one with purely real Ammn(1) and another with purely imaginary Ammn(2) 
% polarizability. The precessing axes are given by real vectors wmmn_i.

an(1,2)=mm12n(f0); an(1,3)=mm13n(f0); an(2,3)=mm23n(f0);
Ammn(1)=sqrt( (real(an(2,3)) )^2 + (real(an(1,3)) )^2 + (real(an(1,2)) )^2 );
Ammn(2)=1i*sqrt( (imag(an(2,3)) )^2 + (imag(an(1,3)) )^2 + (imag(an(1,2)) )^2 );
wmmn(:,1)=1/Ammn(1)*[-real(an(2,3))  real(an(1,3))  -real(an(1,2))]';
wmmn(:,2)=1/imag(Ammn(2))*[-imag(an(2,3))  imag(an(1,3))  -imag(an(1,2))]';

%% Modulus decomposition for alpha_em reciprocal

% the tensor is general and therefore is decomposed to diagonal, symmetric, and antisymmetric parts as: alpha=DemrI + S_emr + A_emr

% diagonal part is realized as 3 equal spiral inclusions oriented along the
% x, y, and z axes. Their amplitudes are equal to Demr
Demr=(em11r(f0)+em22r(f0)+em33r(f0))/3;

% symmetric part is realized as 6 spiral inclusions (with different complex
% aimplitudes Semr_i) oriented along 6 real vectors vemr_i
sr(1,1)=(2*em11r(f0)-em22r(f0)-em33r(f0))/3;
sr(2,2)=(2*em22r(f0)-em11r(f0)-em33r(f0))/3;
sr(3,3)=(2*em33r(f0)-em11r(f0)-em22r(f0))/3;
sr(1,2)=(em12r(f0)+em21r(f0))/2; sr(1,3)=(em13r(f0)+em31r(f0))/2; sr(2,3)=(em23r(f0)+em32r(f0))/2;
vemr(:,1)=[1 0 0]'; vemr(:,2)=[0 1 0]'; vemr(:,3)=[0 0 1]';
for ii=1:2       %this cycle provides the best decomposition (with minimum number of inclusions)
    for jj=1:2
        for kk=1:2
            semr(ii,jj,kk)=abs(sr(1,1)-(-1)^(ii+1)*sr(1,2)-(-1)^(jj+1)*sr(1,3))+abs(sr(2,2)-(-1)^(ii+1)*sr(1,2)...
                -(-1)^(kk+1)*sr(2,3))+abs(sr(3,3)-(-1)^(jj+1)*sr(1,3)-(-1)^(kk+1)*sr(2,3));
        end;
    end
end;
[emmr,Iemr] = min(semr(:));
[iimin, jjmin, kkmin] = ind2sub(size(semr),Iemr);
Semr(1)=sr(1,1)-(-1)^(iimin+1)*sr(1,2)-(-1)^(jjmin+1)*sr(1,3);
Semr(2)=sr(2,2)-(-1)^(iimin+1)*sr(1,2)-(-1)^(kkmin+1)*sr(2,3);
Semr(3)=sr(3,3)-(-1)^(jjmin+1)*sr(1,3)-(-1)^(kkmin+1)*sr(2,3);
Semr(4)=2*(-1)^(iimin+1)*sr(1,2);
Semr(5)=2*(-1)^(jjmin+1)*sr(1,3);
Semr(6)=2*(-1)^(kkmin+1)*sr(2,3);
vemr(:,4)=1/sqrt(2)*[1 (-1)^(iimin+1) 0]'; vemr(:,5)=1/sqrt(2)*[1 0 (-1)^(jjmin+1)]'; vemr(:,6)=1/sqrt(2)*[0 1 (-1)^(kkmin+1)]';

% the antisymmetric part is realized as 2 uniaxial omega inclusions: 
% one with purely real Aemr(1) and another with purely imaginary Aemr(2) 
% polarizability. Their axes are given by real vectors wemr_i.

ar(1,2)=(em12r(f0)-em21r(f0))/2; ar(1,3)=(em13r(f0)-em31r(f0))/2; ar(2,3)=(em23r(f0)-em32r(f0))/2;
Aemr(1)=sqrt( (real(ar(2,3)) )^2 + (real(ar(1,3)) )^2 + (real(ar(1,2)) )^2 );
Aemr(2)=1i*sqrt( (imag(ar(2,3)) )^2 + (imag(ar(1,3)) )^2 + (imag(ar(1,2)) )^2 );
wemr(:,1)=1/Aemr(1)*[-real(ar(2,3))  real(ar(1,3))  -real(ar(1,2))]';
wemr(:,2)=1/imag(Aemr(2))*[-imag(ar(2,3))  imag(ar(1,3))  -imag(ar(1,2))]';

%% Modulus decomposition for alpha_em non-reciprocal

% the tensor is general and therefore is decomposed to diagonal, symmetric, and antisymmetric parts as: alpha=DemnI + S_emn + A_emn

% diagonal part is realized as 3 equal nonreciprocal Tellegen inclusions oriented along the
% x, y, and z axes. Their amplitudes are equal to Demn
Demn=(em11n(f0)+em22n(f0)+em33n(f0))/3;

% symmetric part is realized as 6 nonreciprocal Tellegen inclusions (with different complex
% aimplitudes Semn_i) oriented along 6 real vectors vemn_i
sr(1,1)=(2*em11n(f0)-em22n(f0)-em33n(f0))/3;
sr(2,2)=(2*em22n(f0)-em11n(f0)-em33n(f0))/3;
sr(3,3)=(2*em33n(f0)-em11n(f0)-em22n(f0))/3;
sr(1,2)=(em12n(f0)+em21n(f0))/2; sr(1,3)=(em13n(f0)+em31n(f0))/2; sr(2,3)=(em23n(f0)+em32n(f0))/2;
vemn(:,1)=[1 0 0]'; vemn(:,2)=[0 1 0]'; vemn(:,3)=[0 0 1]';
for ii=1:2       %this cycle provides the best decomposition (with minimum number of inclusions)
    for jj=1:2
        for kk=1:2
            semn(ii,jj,kk)=abs(sr(1,1)-(-1)^(ii+1)*sr(1,2)-(-1)^(jj+1)*sr(1,3))+abs(sr(2,2)-(-1)^(ii+1)*sr(1,2)...
                -(-1)^(kk+1)*sr(2,3))+abs(sr(3,3)-(-1)^(jj+1)*sr(1,3)-(-1)^(kk+1)*sr(2,3));
        end;
    end
end;
[emmn,Iemn] = min(semn(:));
[iimin, jjmin, kkmin] = ind2sub(size(semn),Iemn);
Semn(1)=sr(1,1)-(-1)^(iimin+1)*sr(1,2)-(-1)^(jjmin+1)*sr(1,3);
Semn(2)=sr(2,2)-(-1)^(iimin+1)*sr(1,2)-(-1)^(kkmin+1)*sr(2,3);
Semn(3)=sr(3,3)-(-1)^(jjmin+1)*sr(1,3)-(-1)^(kkmin+1)*sr(2,3);
Semn(4)=2*(-1)^(iimin+1)*sr(1,2);
Semn(5)=2*(-1)^(jjmin+1)*sr(1,3);
Semn(6)=2*(-1)^(kkmin+1)*sr(2,3);
vemn(:,4)=1/sqrt(2)*[1 (-1)^(iimin+1) 0]'; vemn(:,5)=1/sqrt(2)*[1 0 (-1)^(jjmin+1)]'; vemn(:,6)=1/sqrt(2)*[0 1 (-1)^(kkmin+1)]';

% the antisymmetric part is realized as 2 uniaxial "moving" inclusions: 
% one with purely real Aemn(1) and another with purely imaginary Aemn(2) 
% polarizability. Their axes are given by real vectors wemn_i.

ar(1,2)=(em12n(f0)-em21n(f0))/2; ar(1,3)=(em13n(f0)-em31n(f0))/2; ar(2,3)=(em23n(f0)-em32n(f0))/2;
Aemn(1)=sqrt( (real(ar(2,3)) )^2 + (real(ar(1,3)) )^2 + (real(ar(1,2)) )^2 );
Aemn(2)=1i*sqrt( (imag(ar(2,3)) )^2 + (imag(ar(1,3)) )^2 + (imag(ar(1,2)) )^2 );
wemn(:,1)=1/Aemn(1)*[-real(ar(2,3))  real(ar(1,3))  -real(ar(1,2))]';
wemn(:,2)=1/imag(Aemn(2))*[-imag(ar(2,3))  imag(ar(1,3))  -imag(ar(1,2))]';

%% Creating output

theta_ed=180/pi*[acos(veer(3,1)) acos(veer(3,2)) acos(veer(3,3)) acos(veer(3,4)) acos(veer(3,5)) acos(veer(3,6))];
phi_ed=180/pi*[atan2(veer(2,1),veer(1,1)) atan2(veer(2,2),veer(1,2)) atan2(veer(2,3),veer(1,3))...
    atan2(veer(2,4),veer(1,4)) atan2(veer(2,5),veer(1,5)) atan2(veer(2,6),veer(1,6))];
ampl_ed=abs([Deer+Seer(1) Deer+Seer(2) Deer+Seer(3) Seer(4) Seer(5) Seer(6)]);

theta_ped=180/pi*[acos(ween(3,1)) acos(ween(3,2))];
phi_ped=180/pi*[atan2(ween(2,1),ween(1,1)) atan2(ween(2,2),ween(1,2))];
ampl_ped=abs([Aeen(1) Aeen(2)]);

theta_md=180/pi*[acos(vmmr(3,1)) acos(vmmr(3,2)) acos(vmmr(3,3)) acos(vmmr(3,4)) acos(vmmr(3,5)) acos(vmmr(3,6))];
phi_md=180/pi*[atan2(vmmr(2,1),vmmr(1,1)) atan2(vmmr(2,2),vmmr(1,2)) atan2(vmmr(2,3),vmmr(1,3))...
    atan2(vmmr(2,4),vmmr(1,4)) atan2(vmmr(2,5),vmmr(1,5)) atan2(vmmr(2,6),vmmr(1,6))];
ampl_md=abs([Dmmr+Smmr(1) Dmmr+Smmr(2) Dmmr+Smmr(3) Smmr(4) Smmr(5) Smmr(6)]);

theta_pmd=180/pi*[acos(wmmn(3,1)) acos(wmmn(3,2))];
phi_pmd=180/pi*[atan2(wmmn(2,1),wmmn(1,1)) atan2(wmmn(2,2),wmmn(1,2))];
ampl_pmd=abs([Ammn(1) Ammn(2)]);


theta_ch=180/pi*[acos(vemr(3,1)) acos(vemr(3,2)) acos(vemr(3,3)) acos(vemr(3,4)) acos(vemr(3,5)) acos(vemr(3,6))];
phi_ch=180/pi*[atan2(vemr(2,1),vemr(1,1)) atan2(vemr(2,2),vemr(1,2)) atan2(vemr(2,3),vemr(1,3))...
    atan2(vemr(2,4),vemr(1,4)) atan2(vemr(2,5),vemr(1,5)) atan2(vemr(2,6),vemr(1,6))];
ampl_ch=abs([Demr+Semr(1) Demr+Semr(2) Demr+Semr(3) Semr(4) Semr(5) Semr(6)]);
handedness_ch=[sign(-real(Demr+Semr(1))) sign(-real(Demr+Semr(2))) sign(-real(Demr+Semr(3)))...
    sign(-real(Semr(4))) sign(-real(Semr(5))) sign(-real(Semr(6)))];

theta_om=180/pi*[acos(wemr(3,1)) acos(wemr(3,2))];
phi_om=180/pi*[atan2(wemr(2,1),wemr(1,1)) atan2(wemr(2,2),wemr(1,2))];
ampl_om=abs([Aemr(1) Aemr(2)]);


theta_te=180/pi*[acos(vemn(3,1)) acos(vemn(3,2)) acos(vemn(3,3)) acos(vemn(3,4)) acos(vemn(3,5)) acos(vemn(3,6))];
phi_te=180/pi*[atan2(vemn(2,1),vemn(1,1)) atan2(vemn(2,2),vemn(1,2)) atan2(vemn(2,3),vemn(1,3))...
    atan2(vemn(2,4),vemn(1,4)) atan2(vemn(2,5),vemn(1,5)) atan2(vemn(2,6),vemn(1,6))];
ampl_te=abs([Demn+Semn(1) Demn+Semn(2) Demn+Semn(3) Semn(4) Semn(5) Semn(6)]);
orientation_te=[sign(imag(Demn+Semn(1))) sign(imag(Demn+Semn(2))) sign(imag(Demn+Semn(3)))...
    sign(imag(Semn(4))) sign(imag(Semn(5))) sign(imag(Semn(6)))];

theta_mo=180/pi*[acos(wemn(3,1)) acos(wemn(3,2))];
phi_mo=180/pi*[atan2(wemn(2,1),wemn(1,1)) atan2(wemn(2,2),wemn(1,2))];
ampl_mo=abs([Aemn(1) Aemn(2)]);

ampl_max=max([max(ampl_ed) max(ampl_ped) max(ampl_md) max(ampl_pmd) max(ampl_ch) max(ampl_om) max(ampl_te) max(ampl_mo)]);

scale_ed=ampl_ed/ampl_max; scale_ped=ampl_ped/ampl_max;
scale_md=ampl_md/ampl_max; scale_pmd=ampl_pmd/ampl_max;
scale_ch=ampl_ch/ampl_max; scale_om=ampl_om/ampl_max;
scale_te=ampl_te/ampl_max; scale_mo=ampl_mo/ampl_max;

for ss=1:6
    if scale_ed(ss)<1e-5
        scale_ed(ss)=1e-5;
    end
    if scale_md(ss)<1e-5
        scale_md(ss)=1e-5;
    end
    if scale_ch(ss)<1e-5
        scale_ch(ss)=1e-5;
    end
    if scale_te(ss)<1e-5
        scale_te(ss)=1e-5;
    end
end;
for ss=1:2
    if scale_ped(ss)<1e-5
        scale_ped(ss)=1e-5;
    end
    if scale_pmd(ss)<1e-5
        scale_pmd(ss)=1e-5;
    end
    if scale_om(ss)<1e-5
        scale_om(ss)=1e-5;
    end
    if scale_mo(ss)<1e-5
        scale_mo(ss)=1e-5;
    end
end;

%% Finding orientation for maximum chirality and Tellegen effects

max_ch=0; max_te=0;
for tt=0:1:90
    for pp=0:1:360
        g1=sind(tt)*cosd(pp); g2=sind(tt)*sind(pp); g3=cosd(tt); 
        kappa_new=abs(g1^2*(Demr+Semr(1))+g2^2*(Demr+Semr(2))+g3^2*(Demr+Semr(3))+(Semr(4))*(g1*vemr(1,4)+g2*vemr(2,4))^2+...
            (Semr(5))*(g1*vemr(1,5)+g3*vemr(3,5))^2+(Semr(6))*(g2*vemr(2,6)+g3*vemr(3,6))^2);
        if kappa_new > max_ch
           max_ch=kappa_new; theta_opt_ch=tt; phi_opt_ch=pp;
        end
        kappa_total(tt+1,pp+1)=kappa_new;
        
        chi_new=abs(g1^2*(Demn+Semn(1))+g2^2*(Demn+Semn(2))+g3^2*(Demn+Semn(3))+(Semn(4))*(g1*vemn(1,4)+g2*vemn(2,4))^2+...
            (Semn(5))*(g1*vemn(1,5)+g3*vemn(3,5))^2+(Semn(6))*(g2*vemn(2,6)+g3*vemn(3,6))^2);
        if chi_new > max_te
           max_te=chi_new; theta_opt_te=tt; phi_opt_te=pp;
        end
        chi_total(tt+1,pp+1)=chi_new;
    end;
end;
set(0,'DefaultAxesFontName', 'Times New Roman')
set(0,'DefaultAxesFontSize', 14)
[ttt,ppp] = meshgrid(0:1:90,0:1:360);
kappa_total=kappa_total';
chi_total=chi_total';

figure('Name','Chiral polarizability of the meta-atom with respect to different orientations in the spherical coordinate system')      % this figure shows chirality level with respect to the orientation vector g defined by \theta and \phi in a spherical coordinate system
surf(ttt,ppp,kappa_total,'LineStyle',':');
xlabel('$$ \theta \rm{[deg]}$$','Interpreter','latex')
ylabel('$$ \phi \rm{[deg]}$$','Interpreter','latex')
zlabel('$$ |\alpha_{\rm em, rec}^{\rm gg}|  \rm{[m^2 s]}$$','Interpreter','latex')
xlim([0 90]); ylim([0 360]);
set(gca,'XTick',[0:15:90])  
set(gca,'YTick',[0:45:360]) 

figure('Name','Tellegen polarizability of the meta-atom with respect to different orientations in the spherical coordinate system')      % this figure shows Tellegen coupling level with respect to the orientation vector g defined by \theta and \phi in a spherical coordinate system
surf(ttt,ppp,chi_total,'LineStyle',':');
xlabel('$$ \theta \rm{[deg]}$$','Interpreter','latex')
ylabel('$$ \phi \rm{[deg]}$$','Interpreter','latex')
zlabel('$$ |\alpha_{\rm em, nonr}^{\rm gg}|  \rm{[m^2 s]}$$','Interpreter','latex')
xlim([0 90]); ylim([0 360]);
set(gca,'XTick',[0:15:90])  
set(gca,'YTick',[0:45:360])


%  theta_opt_ch=0; phi_opt_ch=0;
g1=sind(theta_opt_ch)*cosd(phi_opt_ch); g2=sind(theta_opt_ch)*sind(phi_opt_ch); g3=cosd(theta_opt_ch); 
max_ch_complex=g1^2*(Demr+Semr(1))+g2^2*(Demr+Semr(2))+g3^2*(Demr+Semr(3))+(Semr(4))*(g1*vemr(1,4)+g2*vemr(2,4))^2+...
            (Semr(5))*(g1*vemr(1,5)+g3*vemr(3,5))^2+(Semr(6))*(g2*vemr(2,6)+g3*vemr(3,6))^2;
handedness_ch_max=sign(-real(max_ch_complex));
scale_ch_max=abs(max_ch_complex)/ampl_max;


g1=sind(theta_opt_te)*cosd(phi_opt_te); g2=sind(theta_opt_te)*sind(phi_opt_te); g3=cosd(theta_opt_te); 
max_te_complex=g1^2*(Demn+Semn(1))+g2^2*(Demn+Semn(2))+g3^2*(Demn+Semn(3))+(Semn(4))*(g1*vemn(1,4)+g2*vemn(2,4))^2+...
            (Semn(5))*(g1*vemn(1,5)+g3*vemn(3,5))^2+(Semn(6))*(g2*vemn(2,6)+g3*vemn(3,6))^2;
orientation_te_max=sign(imag(max_te_complex));
scale_te_max=abs(max_te_complex)/ampl_max;

if scale_ch_max<1e-5
        scale_ch_max=1e-5;
end
if scale_te_max<1e-5
        scale_te_max=1e-5;
end

%% Generating script for ANSYS HFSS project "Materiatronics_modules.aedt"

fileID = fopen([name '.py'],'w');
ProjectName='Materiatronics_modules';
fprintf(fileID, 'import ScriptEnv\n');
fprintf(fileID, 'ScriptEnv.Initialize("Ansoft.ElectronicsDesktop")\n');
fprintf(fileID, 'oDesktop.RestoreWindow()\n');
fprintf(fileID, 'oProject = oDesktop.SetActiveProject("%s")\n',ProjectName);
fprintf(fileID, 'oDesign = oProject.SetActiveDesign("HFSSDesign1")\n');
fprintf(fileID, 'oEditor = oDesign.SetActiveEditor("3D Modeler")\n');


% ee reciprocal
for ic=1:6
fprintf(fileID, 'oEditor.ChangeProperty(["NAME:AllTabs",["NAME:Geometry3DCmdTab",["NAME:PropServers","Electric_dipole_%s:Scale:1"],["NAME:ChangedProps",',num2str(ic));
fprintf(fileID, '["NAME:Scale of X","Value:=", "%s"],',num2str(scale_ed(ic)));
fprintf(fileID, '["NAME:Scale of Y","Value:=", "%s"],',num2str(scale_ed(ic)));
fprintf(fileID, '["NAME:Scale of Z","Value:=", "%s"]]]])\n',num2str(scale_ed(ic)));
fprintf(fileID, 'oEditor.ChangeProperty(["NAME:AllTabs",["NAME:Geometry3DCmdTab",["NAME:PropServers","Electric_dipole_%s:Rotate:1"],',num2str(ic));
fprintf(fileID, '["NAME:ChangedProps",["NAME:Angle",	"Value:="," %sdeg"]]]])\n',num2str(theta_ed(ic)));
fprintf(fileID, 'oEditor.ChangeProperty(["NAME:AllTabs",["NAME:Geometry3DCmdTab",["NAME:PropServers","Electric_dipole_%s:Rotate:2"],',num2str(ic));
fprintf(fileID, '["NAME:ChangedProps",["NAME:Angle",	"Value:="," %sdeg"]]]])\n',num2str(phi_ed(ic)));
end;

% ee non-reciprocal
fprintf(fileID, 'oEditor.ChangeProperty(["NAME:AllTabs",["NAME:Geometry3DCmdTab",["NAME:PropServers","Precessing_electric_dipole_real:Scale:2"],["NAME:ChangedProps",');
fprintf(fileID, '["NAME:Scale of X","Value:=", "%s"],',num2str(scale_ped(1)));
fprintf(fileID, '["NAME:Scale of Y","Value:=", "%s"],',num2str(scale_ped(1)));
fprintf(fileID, '["NAME:Scale of Z","Value:=", "%s"]]]])\n',num2str(scale_ped(1)));
fprintf(fileID, 'oEditor.ChangeProperty(["NAME:AllTabs",["NAME:Geometry3DCmdTab",["NAME:PropServers","Precessing_electric_dipole_real:Rotate:1"],');
fprintf(fileID, '["NAME:ChangedProps",["NAME:Angle",	"Value:="," %sdeg"]]]])\n',num2str(round(theta_ped(1))));
fprintf(fileID, 'oEditor.ChangeProperty(["NAME:AllTabs",["NAME:Geometry3DCmdTab",["NAME:PropServers","Precessing_electric_dipole_real:Rotate:2"],');
fprintf(fileID, '["NAME:ChangedProps",["NAME:Angle",	"Value:="," %sdeg"]]]])\n',num2str(round(phi_ped(1))));

fprintf(fileID, 'oEditor.ChangeProperty(["NAME:AllTabs",["NAME:Geometry3DCmdTab",["NAME:PropServers","Precessing_electric_dipole_imaginary:Scale:2"],["NAME:ChangedProps",');
fprintf(fileID, '["NAME:Scale of X","Value:=", "%s"],',num2str(scale_ped(2)));
fprintf(fileID, '["NAME:Scale of Y","Value:=", "%s"],',num2str(scale_ped(2)));
fprintf(fileID, '["NAME:Scale of Z","Value:=", "%s"]]]])\n',num2str(scale_ped(2)));
fprintf(fileID, 'oEditor.ChangeProperty(["NAME:AllTabs",["NAME:Geometry3DCmdTab",["NAME:PropServers","Precessing_electric_dipole_imaginary:Rotate:1"],');
fprintf(fileID, '["NAME:ChangedProps",["NAME:Angle",	"Value:="," %sdeg"]]]])\n',num2str(round(theta_ped(2))));
fprintf(fileID, 'oEditor.ChangeProperty(["NAME:AllTabs",["NAME:Geometry3DCmdTab",["NAME:PropServers","Precessing_electric_dipole_imaginary:Rotate:2"],');
fprintf(fileID, '["NAME:ChangedProps",["NAME:Angle",	"Value:="," %sdeg"]]]])\n',num2str(round(phi_ped(2))));

% mm reciprocal
for ic=1:6
fprintf(fileID, 'oEditor.ChangeProperty(["NAME:AllTabs",["NAME:Geometry3DCmdTab",["NAME:PropServers","Cylinder%s:Scale:1"],["NAME:ChangedProps",',num2str(ic*2+2));
fprintf(fileID, '["NAME:Scale of X","Value:=", "%s"],',num2str(scale_md(ic)));
fprintf(fileID, '["NAME:Scale of Y","Value:=", "%s"],',num2str(scale_md(ic)));
fprintf(fileID, '["NAME:Scale of Z","Value:=", "%s"]]]])\n',num2str(scale_md(ic)));
fprintf(fileID, 'oEditor.ChangeProperty(["NAME:AllTabs",["NAME:Geometry3DCmdTab",["NAME:PropServers","Cylinder%s:Rotate:1"],',num2str(ic*2+2));
fprintf(fileID, '["NAME:ChangedProps",["NAME:Angle",	"Value:="," %sdeg"]]]])\n',num2str(theta_md(ic)));
fprintf(fileID, 'oEditor.ChangeProperty(["NAME:AllTabs",["NAME:Geometry3DCmdTab",["NAME:PropServers","Cylinder%s:Rotate:2"],',num2str(ic*2+2));
fprintf(fileID, '["NAME:ChangedProps",["NAME:Angle",	"Value:="," %sdeg"]]]])\n',num2str(phi_md(ic)));
end;
 
% mm non-reciprocal
fprintf(fileID, 'oEditor.ChangeProperty(["NAME:AllTabs",["NAME:Geometry3DCmdTab",["NAME:PropServers","Precessing_magnetic_dipole_real:Scale:2"],["NAME:ChangedProps",');
fprintf(fileID, '["NAME:Scale of X","Value:=", "%s"],',num2str(scale_pmd(1)));
fprintf(fileID, '["NAME:Scale of Y","Value:=", "%s"],',num2str(scale_pmd(1)));
fprintf(fileID, '["NAME:Scale of Z","Value:=", "%s"]]]])\n',num2str(scale_pmd(1)));
fprintf(fileID, 'oEditor.ChangeProperty(["NAME:AllTabs",["NAME:Geometry3DCmdTab",["NAME:PropServers","Precessing_magnetic_dipole_real:Rotate:1"],');
fprintf(fileID, '["NAME:ChangedProps",["NAME:Angle",	"Value:="," %sdeg"]]]])\n',num2str(round(theta_pmd(1))));
fprintf(fileID, 'oEditor.ChangeProperty(["NAME:AllTabs",["NAME:Geometry3DCmdTab",["NAME:PropServers","Precessing_magnetic_dipole_real:Rotate:2"],');
fprintf(fileID, '["NAME:ChangedProps",["NAME:Angle",	"Value:="," %sdeg"]]]])\n',num2str(round(phi_pmd(1))));

fprintf(fileID, 'oEditor.ChangeProperty(["NAME:AllTabs",["NAME:Geometry3DCmdTab",["NAME:PropServers","Precessing_magnetic_dipole_imaginary:Scale:2"],["NAME:ChangedProps",');
fprintf(fileID, '["NAME:Scale of X","Value:=", "%s"],',num2str(scale_pmd(2)));
fprintf(fileID, '["NAME:Scale of Y","Value:=", "%s"],',num2str(scale_pmd(2)));
fprintf(fileID, '["NAME:Scale of Z","Value:=", "%s"]]]])\n',num2str(scale_pmd(2)));
fprintf(fileID, 'oEditor.ChangeProperty(["NAME:AllTabs",["NAME:Geometry3DCmdTab",["NAME:PropServers","Precessing_magnetic_dipole_imaginary:Rotate:1"],');
fprintf(fileID, '["NAME:ChangedProps",["NAME:Angle",	"Value:="," %sdeg"]]]])\n',num2str(round(theta_pmd(2))));
fprintf(fileID, 'oEditor.ChangeProperty(["NAME:AllTabs",["NAME:Geometry3DCmdTab",["NAME:PropServers","Precessing_magnetic_dipole_imaginary:Rotate:2"],');
fprintf(fileID, '["NAME:ChangedProps",["NAME:Angle",	"Value:="," %sdeg"]]]])\n',num2str(round(phi_pmd(2))));

% em reciprocal
for ic=1:6
fprintf(fileID, 'oEditor.ChangeProperty(["NAME:AllTabs",["NAME:Geometry3DCmdTab",["NAME:PropServers","EquationCurve%s:Scale:1"],["NAME:ChangedProps",',num2str(ic));
fprintf(fileID, '["NAME:Scale of X","Value:=", "%s"],',num2str(scale_ch(ic)));
fprintf(fileID, '["NAME:Scale of Y","Value:=", "%s"],',num2str(scale_ch(ic)));
fprintf(fileID, '["NAME:Scale of Z","Value:=", "%s"]]]])\n',num2str(scale_ch(ic)));
fprintf(fileID, 'oEditor.ChangeProperty(["NAME:AllTabs",["NAME:Geometry3DCmdTab",["NAME:PropServers","EquationCurve%s:Rotate:1"],',num2str(ic));
fprintf(fileID, '["NAME:ChangedProps",["NAME:Angle",	"Value:="," %sdeg"]]]])\n',num2str(theta_ch(ic)));
fprintf(fileID, 'oEditor.ChangeProperty(["NAME:AllTabs",["NAME:Geometry3DCmdTab",["NAME:PropServers","EquationCurve%s:Rotate:2"],',num2str(ic));
fprintf(fileID, '["NAME:ChangedProps",["NAME:Angle",	"Value:="," %sdeg"]]]])\n',num2str(phi_ch(ic)));
fprintf(fileID, 'oEditor.ChangeProperty(["NAME:AllTabs",["NAME:Geometry3DCmdTab",["NAME:PropServers","EquationCurve%s:CreateEquationCurve:1"],',num2str(ic));
fprintf(fileID, '["NAME:ChangedProps",["NAME:Y(_t)","Value:=","%slambda/25*sin(_t)"]]]])\n',[num2str(handedness_ch(ic)) '*']);
end;

for ic=1:2
fprintf(fileID, 'oEditor.ChangeProperty(["NAME:AllTabs",["NAME:Geometry3DCmdTab",["NAME:PropServers","Cylinder%s:Scale:1"],["NAME:ChangedProps",',num2str(ic*2+22));
fprintf(fileID, '["NAME:Scale of X","Value:=", "%s"],',num2str(scale_om(ic)));
fprintf(fileID, '["NAME:Scale of Y","Value:=", "%s"],',num2str(scale_om(ic)));
fprintf(fileID, '["NAME:Scale of Z","Value:=", "%s"]]]])\n',num2str(scale_om(ic)));
fprintf(fileID, 'oEditor.ChangeProperty(["NAME:AllTabs",["NAME:Geometry3DCmdTab",["NAME:PropServers","Cylinder%s:Rotate:3"],',num2str(ic*2+22));
fprintf(fileID, '["NAME:ChangedProps",["NAME:Angle",	"Value:="," %sdeg"]]]])\n',num2str(round(theta_om(ic))));
fprintf(fileID, 'oEditor.ChangeProperty(["NAME:AllTabs",["NAME:Geometry3DCmdTab",["NAME:PropServers","Cylinder%s:Rotate:4"],',num2str(ic*2+22));
fprintf(fileID, '["NAME:ChangedProps",["NAME:Angle",	"Value:="," %sdeg"]]]])\n',num2str(round(phi_om(ic))));
end;


% em non-reciprocal
for ic=1:6
fprintf(fileID, 'oEditor.ChangeProperty(["NAME:AllTabs",["NAME:Geometry3DCmdTab",["NAME:PropServers","striptel%s:Scale:1"],["NAME:ChangedProps",',num2str(ic));
fprintf(fileID, '["NAME:Scale of X","Value:=", "%s"],',num2str(scale_te(ic)));
fprintf(fileID, '["NAME:Scale of Y","Value:=", "%s"],',num2str(scale_te(ic)));
fprintf(fileID, '["NAME:Scale of Z","Value:=", "%s"]]]])\n',num2str(scale_te(ic)));
fprintf(fileID, 'oEditor.ChangeProperty(["NAME:AllTabs",["NAME:Geometry3DCmdTab",["NAME:PropServers","striptel%s:Rotate:1"],',num2str(ic));
fprintf(fileID, '["NAME:ChangedProps",["NAME:Angle",	"Value:="," %sdeg"]]]])\n',num2str(theta_te(ic)));
fprintf(fileID, 'oEditor.ChangeProperty(["NAME:AllTabs",["NAME:Geometry3DCmdTab",["NAME:PropServers","striptel%s:Rotate:2"],',num2str(ic));
fprintf(fileID, '["NAME:ChangedProps",["NAME:Angle",	"Value:="," %sdeg"]]]])\n',num2str(phi_te(ic)));

fprintf(fileID, 'oEditor.ChangeProperty(["NAME:AllTabs",["NAME:Geometry3DCmdTab",["NAME:PropServers","striptel%s:CreateCylinder:1"],["NAME:ChangedProps",',num2str(ic));
fprintf(fileID, '["NAME:Center Position","X:=","0","Y:=","%slambda/25","Z:=", "-lambda/5/2"]]]])\n',[num2str(orientation_te(ic)) '*']);
end;

fprintf(fileID, 'oEditor.ChangeProperty(["NAME:AllTabs",["NAME:Geometry3DCmdTab",["NAME:PropServers","Moving_real:Scale:2"],["NAME:ChangedProps",');
fprintf(fileID, '["NAME:Scale of X","Value:=", "%s"],',num2str(scale_mo(1)));
fprintf(fileID, '["NAME:Scale of Y","Value:=", "%s"],',num2str(scale_mo(1)));
fprintf(fileID, '["NAME:Scale of Z","Value:=", "%s"]]]])\n',num2str(scale_mo(1)));
fprintf(fileID, 'oEditor.ChangeProperty(["NAME:AllTabs",["NAME:Geometry3DCmdTab",["NAME:PropServers","Moving_real:Rotate:1"],');
fprintf(fileID, '["NAME:ChangedProps",["NAME:Angle",	"Value:="," %sdeg"]]]])\n',num2str(round(theta_mo(1))));
fprintf(fileID, 'oEditor.ChangeProperty(["NAME:AllTabs",["NAME:Geometry3DCmdTab",["NAME:PropServers","Moving_real:Rotate:2"],');
fprintf(fileID, '["NAME:ChangedProps",["NAME:Angle",	"Value:="," %sdeg"]]]])\n',num2str(round(phi_mo(1))));

fprintf(fileID, 'oEditor.ChangeProperty(["NAME:AllTabs",["NAME:Geometry3DCmdTab",["NAME:PropServers","Moving_imag:Scale:2"],["NAME:ChangedProps",');
fprintf(fileID, '["NAME:Scale of X","Value:=", "%s"],',num2str(scale_mo(2)));
fprintf(fileID, '["NAME:Scale of Y","Value:=", "%s"],',num2str(scale_mo(2)));
fprintf(fileID, '["NAME:Scale of Z","Value:=", "%s"]]]])\n',num2str(scale_mo(2)));
fprintf(fileID, 'oEditor.ChangeProperty(["NAME:AllTabs",["NAME:Geometry3DCmdTab",["NAME:PropServers","Moving_imag:Rotate:1"],');
fprintf(fileID, '["NAME:ChangedProps",["NAME:Angle",	"Value:="," %sdeg"]]]])\n',num2str(round(theta_mo(2))));
fprintf(fileID, 'oEditor.ChangeProperty(["NAME:AllTabs",["NAME:Geometry3DCmdTab",["NAME:PropServers","Moving_imag:Rotate:2"],');
fprintf(fileID, '["NAME:ChangedProps",["NAME:Angle",	"Value:="," %sdeg"]]]])\n',num2str(round(phi_mo(2))));


% orientation for maximum chirality and Tellegen effects
fprintf(fileID, 'oEditor.ChangeProperty(["NAME:AllTabs",["NAME:Geometry3DCmdTab",["NAME:PropServers","EquationCurve7:Scale:1"],["NAME:ChangedProps",');
fprintf(fileID, '["NAME:Scale of X","Value:=", "%s"],',num2str(scale_ch_max));
fprintf(fileID, '["NAME:Scale of Y","Value:=", "%s"],',num2str(scale_ch_max));
fprintf(fileID, '["NAME:Scale of Z","Value:=", "%s"]]]])\n',num2str(scale_ch_max));
fprintf(fileID, 'oEditor.ChangeProperty(["NAME:AllTabs",["NAME:Geometry3DCmdTab",["NAME:PropServers","EquationCurve7:Rotate:1"],');
fprintf(fileID, '["NAME:ChangedProps",["NAME:Angle",	"Value:="," %sdeg"]]]])\n',num2str(theta_opt_ch));
fprintf(fileID, 'oEditor.ChangeProperty(["NAME:AllTabs",["NAME:Geometry3DCmdTab",["NAME:PropServers","EquationCurve7:Rotate:2"],');
fprintf(fileID, '["NAME:ChangedProps",["NAME:Angle",	"Value:="," %sdeg"]]]])\n',num2str(phi_opt_ch));
fprintf(fileID, 'oEditor.ChangeProperty(["NAME:AllTabs",["NAME:Geometry3DCmdTab",["NAME:PropServers","EquationCurve7:CreateEquationCurve:1"],');
fprintf(fileID, '["NAME:ChangedProps",["NAME:Y(_t)","Value:=","%slambda/25*sin(_t)"]]]])\n',[num2str(handedness_ch_max) '*']);
 
fprintf(fileID, 'oEditor.ChangeProperty(["NAME:AllTabs",["NAME:Geometry3DCmdTab",["NAME:PropServers","striptelmax:Scale:1"],["NAME:ChangedProps",');
fprintf(fileID, '["NAME:Scale of X","Value:=", "%s"],',num2str(scale_te_max));
fprintf(fileID, '["NAME:Scale of Y","Value:=", "%s"],',num2str(scale_te_max));
fprintf(fileID, '["NAME:Scale of Z","Value:=", "%s"]]]])\n',num2str(scale_te_max));
fprintf(fileID, 'oEditor.ChangeProperty(["NAME:AllTabs",["NAME:Geometry3DCmdTab",["NAME:PropServers","striptelmax:Rotate:1"],');
fprintf(fileID, '["NAME:ChangedProps",["NAME:Angle",	"Value:="," %sdeg"]]]])\n',num2str(theta_opt_te));
fprintf(fileID, 'oEditor.ChangeProperty(["NAME:AllTabs",["NAME:Geometry3DCmdTab",["NAME:PropServers","striptelmax:Rotate:2"],');
fprintf(fileID, '["NAME:ChangedProps",["NAME:Angle",	"Value:="," %sdeg"]]]])\n',num2str(phi_opt_te));
fprintf(fileID, 'oEditor.ChangeProperty(["NAME:AllTabs",["NAME:Geometry3DCmdTab",["NAME:PropServers","striptelmax:CreateCylinder:1"],["NAME:ChangedProps",');
fprintf(fileID, '["NAME:Center Position","X:=","0","Y:=","%slambda/25","Z:=", "-lambda/5/2"]]]])\n',[num2str(orientation_te_max) '*']);

fclose(fileID);

%% Generating a data table

row_names_symm = {'1';'2';'3';'4';'5';'6'};
row_names_antisymm = {'1';'2'};
col_names_symm = {'s_eer';'theta_eer';'phi_eer';'s_mmr';'theta_mmr';'phi_mmr';'kappa';'theta_kappa';'phi_kappa';'chi';'theta_chi';'phi_chi';'H_kappa';'H_chi'};
col_names_antisymm = {'a_een';'theta_een';'phi_een';'a_mmn';'theta_mmn';'phi_mmn';'Omega';'theta_Omega';'phi_Omega';'v';'theta_v';'phi_v'};

strdes='This table contains all the data of the modular decomposition of your meta-atom. The following definitions are used in the table:';
streer='$$\overline{\overline{\alpha}}_{\rm eer}=\sum_{i=1}^6 s_{{\rm eer,}i}   \overline{\nu}_{{\rm eer,}i} \overline{\nu}_{{\rm eer,}i}  $$';
strmmr='$$\overline{\overline{\alpha}}_{\rm mmr}=\sum_{i=1}^6 s_{{\rm mmr,}i}   \overline{\nu}_{{\rm mmr,}i} \overline{\nu}_{{\rm mmr,}i}  $$';
streen='$$\overline{\overline{\alpha}}_{\rm een}= a_{\rm een, 1}  \overline{w}_{\rm een, 1}\times \overline{\overline{I}} + j a_{\rm een, 2}  \overline{w}_{\rm een, 2}\times \overline{\overline{I}}$$';
strmmn='$$\overline{\overline{\alpha}}_{\rm mmn}= a_{\rm mmn, 1}  \overline{w}_{\rm mmn, 1}\times \overline{\overline{I}} + j a_{\rm mmn, 2}  \overline{w}_{\rm mmn, 2}\times \overline{\overline{I}}$$';
stremr='$$\overline{\overline{\alpha}}_{\rm emr}= \sum_{i=1}^6 \kappa_i   \overline{\nu}_{\kappa ,i} \overline{\nu}_{\kappa ,i}+ \Omega_1  \overline{w}_{\Omega  1}\times \overline{\overline{I}} + j \Omega_2  \overline{w}_{\Omega 2}\times \overline{\overline{I}}$$';
stremn='$$\overline{\overline{\alpha}}_{\rm emn}= \sum_{i=1}^6 \chi_i   \overline{\nu}_{\chi ,i} \overline{\nu}_{\chi ,i}+ v_1  \overline{w}_{v  1}\times \overline{\overline{I}} + j v_2  \overline{w}_{v 2}\times \overline{\overline{I}}$$';
strdes2=['All the above unit vectors are defined by the corresponding' newline 'azimuth and polar angles. For example:' newline];
strvectors='$$ \overline{\nu}_{{\rm eer,}i}=\sin\theta_{{\rm eer},i} \cos\phi_{{\rm ee},i} \overline{x}+ \sin\theta_{{\rm ee},i} \sin\phi_{{\rm ee},i} \overline{y}+\cos\theta_{{\rm ee},i}  \overline{z}.$$';

anot1=[strdes,' ', streer,'$$\hspace{5.2cm}$$',strmmr,' ',streen,'$$\hspace{2.8cm}$$',strmmn];
anot2=[stremr newline stremn];
anot3=[strdes2,' ',strvectors,newline,'In the table below, $$H_{\kappa}$$ and $$H_{\chi}$$ denote handedness of chiral' newline 'and Tellegen modules, respectively.' newline 'All amplitude coefficients ($$s_{{\rm eer,}i}$$, $$a_{{\rm eer,}i}$$, ...) were normalized' newline 'by value ' num2str(ampl_max) '. Only the absolute values' newline 'of the normalized coefficients are    shown in the table.'];

for hh=1:6
if handedness_ch(hh)==1
       strhand_ch(hh)="Right";
elseif handedness_ch(hh)==-1 
       strhand_ch(hh)="Left";
end;
if orientation_te(hh)==1
       strhand_te(hh)="Right";
elseif orientation_te(hh)==-1 
       strhand_te(hh)="Left";
end;
end;
symm={scale_ed(1) theta_ed(1) phi_ed(1) scale_md(1) theta_md(1) phi_md(1) scale_ch(1) theta_ch(1) phi_ch(1) scale_te(1) theta_te(1) phi_te(1) convertStringsToChars(strhand_ch(1)) convertStringsToChars(strhand_te(1));
    scale_ed(2) theta_ed(2) phi_ed(2) scale_md(2) theta_md(2) phi_md(2) scale_ch(2) theta_ch(2) phi_ch(2) scale_te(2) theta_te(2) phi_te(2) convertStringsToChars(strhand_ch(2)) convertStringsToChars(strhand_te(2));
    scale_ed(3) theta_ed(3) phi_ed(3) scale_md(3) theta_md(3) phi_md(3) scale_ch(3) theta_ch(3) phi_ch(3) scale_te(3) theta_te(3) phi_te(3) convertStringsToChars(strhand_ch(3)) convertStringsToChars(strhand_te(3));
    scale_ed(4) theta_ed(4) phi_ed(4) scale_md(4) theta_md(4) phi_md(4) scale_ch(4) theta_ch(4) phi_ch(4) scale_te(4) theta_te(4) phi_te(4) convertStringsToChars(strhand_ch(4)) convertStringsToChars(strhand_te(4));
    scale_ed(5) theta_ed(5) phi_ed(5) scale_md(5) theta_md(5) phi_md(5) scale_ch(5) theta_ch(5) phi_ch(5) scale_te(5) theta_te(5) phi_te(5) convertStringsToChars(strhand_ch(5)) convertStringsToChars(strhand_te(5));
    scale_ed(6) theta_ed(6) phi_ed(6) scale_md(6) theta_md(6) phi_md(6) scale_ch(6) theta_ch(6) phi_ch(6) scale_te(6) theta_te(6) phi_te(6) convertStringsToChars(strhand_ch(6)) convertStringsToChars(strhand_te(6)) };
antisymm=[scale_ped(1) theta_ped(1) phi_ped(1) scale_pmd(1) theta_pmd(1) phi_pmd(1) scale_om(1) theta_om(1) phi_om(1) scale_mo(1) theta_mo(1) phi_mo(1);
    scale_ped(2) theta_ped(2) phi_ped(2) scale_pmd(2) theta_pmd(2) phi_pmd(2) scale_om(2) theta_om(2) phi_om(2) scale_mo(2) theta_mo(2) phi_mo(2) ];

figure1 = figure('Name','Table','units','normalized','outerposition',[0.1 0.1 0.8 0.8]);
annotation(figure1,'textbox',[0.1 0.45 0.5 0.5],'String',anot1,'Interpreter','latex','FitBoxToText','off','EdgeColor','none','FontSize',15);
annotation(figure1,'textbox',[0.1 0.25 0.5 0.5],'String',anot2,'Interpreter','latex','FitBoxToText','off','EdgeColor','none','FontSize',15);
annotation(figure1,'textbox',[0.43 0.23 0.5 0.5],'String',anot3,'Interpreter','latex','FitBoxToText','on','EdgeColor','k','FontSize',15);
cs=70;
col_width1={cs cs cs cs cs cs cs cs cs cs cs cs cs cs}; col_width2={cs cs cs cs cs cs cs cs cs cs cs cs};
unit1=uitable('ColumnWidth',col_width1,'Data', symm, 'RowName',row_names_symm,'ColumnName', col_names_symm, 'Position', [130 130 1020 130]);
unit2=uitable('ColumnWidth',col_width2,'Data', antisymm, 'RowName',row_names_antisymm,'ColumnName', col_names_antisymm, 'Position', [130 60 880 60]);
