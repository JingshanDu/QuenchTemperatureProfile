% Quench Temperature Profile Estimation
% Jingshan S. Du | Northwestern University | du@u.northwestern.edu

% A script that estimates the heat transfer / cooling profile of a solid
% wafer/plate in a continuous binary mixed gas flow.

% In this example, cooling of a silicon or glassy carbon plate in a
% hydrogen:argon mixture flow in a tube is calculated.

% The resulting plot shows two curves: (1) plate temperature (assuming this
% is uniform across the plate), and (2) surface flow temperature (defined
% as the average of the plate and gas temperature).

% Notes:
% 1. Heat transfer rates through different mechanisms are added up
% 2. Forced convection has a laminar flow
% 3. Data sources and calculation methods are referenced in each function
% 4. Solid conduction has not been implemented in this script

% Directions:
% Customize the gas flow and physical properties of the plates/gas to
% satisfy your system conditions. Enable or disable different heat transfer
% mechanisms in the Parameter Setup section. Then this code should run and
% plot the temperature profile.

clear all;
global R g Tm D Qm um Tp w l d Lc flagPlateMaterial MSi MC rhoSi rhoC ArToH2VolRatio xH2 xAr MH2 MAr TbH2 TbAr
%% === Parameter Setup ===
flagForcedConvectionEnabled = 1;
flagNaturalConvectionEnabled = 1;
flagRadiationEnabled = 1;
flagPlateMaterial = 1; % 1 = silicon, 2 = glassy carbon
flagTestMode = 0; % Enable to only calculate q_fc, q_nc, and q_r
flagVerbose = 0; % Enable to write each timestep to the console
%% = Useful Tools and Constants =
R = 8.314; % J/K`mol
sigma = 5.67E-8; % W/m2`K4
g = 9.80665; % m/s2

%% = Gas Flow Setup =
Tm = TempCtoK(20); % input gas temperature
if flagPlateMaterial == 1
    D = 0.02; % tube inner diameter, m
    QH2sccm = 100; % sccm = cm3/min
    QArsccm = 200;
elseif flagPlateMaterial == 2
    D = 0.027; % tube inner diameter, m
    QH2sccm = 100; % sccm = cm3/min
    QArsccm = 2000;
end
Qmsccm = QH2sccm + QArsccm;
ArToH2VolRatio = QArsccm/QH2sccm; % Volumetric ratio of flow rate
Qm = Qmsccm/1E6/60; % m3/s
um = Qm/pi/(D/2)^2; % flow velocity, m/s

%% = Initial Plate Condition Setup =
Tp = TempCtoK(600); % initial plate temperature
vTf = (Tp+Tm)/2; % initial surface flow temperature

%% = Physical Properties of the Plate =
if flagPlateMaterial == 1
    % silicon dimensions
    w = 0.003; % m
    l = 0.003; % m
    d = 0.0002; % m
    e = 0.95; % emissivity
elseif flagPlateMaterial == 2
    % glassy carbon dimensions
    w = 0.025;
    l = 0.025;
    d = 0.001;
    e = 0.85; % emissivity
end
Lc = w*l/2/(w+l);
MSi = 28.0855;
MC = 12.0107;
rhoSi = 2328; % kg/m3
rhoC = 1.42E3; % from HTW, kg/m3

%% = Physical Properties of Binary Gas Mixture =
xH2 = 1/(1+ArToH2VolRatio);
xAr = ArToH2VolRatio/(1+ArToH2VolRatio);
MH2 = 2.0158; % Molar mass g/mol
MAr = 39.948; % g/mol
TbH2 = TempCtoK(-252.9); % boiling point
TbAr = TempCtoK(-185.8);

%% === Main Entry ===
% Set up numerical steps:
dtime = 1E-2; % s
% Set up stop condition:
limitdT = 10; % stops if Tp is within limitdT K above Tm.
% init a few vars
i = 1; % loop index
dTp = 0;
if flagPlateMaterial == 1
    rhop = rhoSi;
elseif flagPlateMaterial == 2
    rhop = rhoC;
end
time = 0;
h_fc = 0;
h_nc = 0;
h_c = 0;
q_c = 0;
q_r = 0;
% result matrix:
result = [time,Tp,vTf];
% loop:
if flagTestMode
    q_fc = Nu_fc(Tp)*km(Tf(Tp))*w*(Tp-Tm)
    q_nc = Nu_fc(Tp)*km(Tf(Tp))/Lc*w*l*(Tp-Tm)
    q_r = sigma*w*l*e*(Tp^4-Tm^4)
else
    while Tp > Tm + limitdT
        i = i + 1;
        if flagForcedConvectionEnabled
            h_fc = Nu_fc(Tp)*km(Tf(Tp))/l;
            if flagPlateMaterial == 2
                h_fc = h_fc * 2; % two surfaces for GC
            end
        end
        if flagNaturalConvectionEnabled
            h_nc = Nu_fc(Tp)*km(Tf(Tp))/Lc;
        end
        if flagRadiationEnabled
            q_r = sigma*w*l*e*(Tp^4-Tm^4); 
            if flagPlateMaterial == 2
                q_r = q_r * 2; % two surfaces for GC
            end
        end
        h_c = (h_fc^3.2+h_nc^3.2)^(1/3.2);
        q_c = h_c*w*l*(Tp-Tm);
        dTp = (q_c+q_r)/(Cpp(Tp)*rhop*w*l*d)*dtime;
        Tp = Tp - dTp;
        vTf = (Tm + Tp)/2;
        time = (i-1)*dtime;
        result(i,:) = [time,Tp,vTf];
        if flagVerbose
            fprintf('Step %05d, time %02.4f s, Tp %03.2f K.\n',i,time,Tp)
        end
    end
    % plot:
    plot(result(:,1),result(:,2),result(:,1),result(:,3))
    xlabel('Time (s)');
    ylabel('Plate temperature (K)');
end



%% === Functions ===
%% = func: Useful Tools =
function TK = TempCtoK(TC)
TK = TC + 273.15;
end

%% = func: Physical Properties of the Plate =
function Cp = Cpp(T) % J/kg`K
global MSi MC flagPlateMaterial
if flagPlateMaterial == 1
    % DOI:10.18434/T42S31
    t = T/1000;
    A = 22.81719;
    B = 3.899510;
    C = -0.082885;
    D = 0.042111;
    E = -0.354063;
    Cp = A + B*t + C*t^2 + D*t^3 + E/t^2; % J/mol`K
    Cp = Cp / MSi * 1000;
elseif flagPlateMaterial == 2
    % 2nd poly fit of DOI:10.7209/tanso.1971.44
    Cp = 1.8851 + 0.00822*T - 5.20848E-6*T^2; % cal/mol`K
    Cp = 4.184 *Cp / MC * 1000;
end
end

%% = func: Physical Properties of Binary Gas Mixture =
function rho = rhoH2(T) % kg/m3
rho = 0.0893E-3 / T * TempCtoK(25); % assuming ideal gas
end

function rho = rhoAr(T) % kg/m3
rho = 1.784 / T * TempCtoK(25); % assuming ideal gas
end

function rho = rhom(T) % kg/m3
global xH2 xAr
rho = (rhoH2(T)*xH2 + rhoAr(T)*xAr);
end

function Cp = CpH2(T) % J/kg`K
global MH2
% DOI:10.18434/T42S31
t = T/1000;
A = 33.066178;
B = -11.363417;
C = 11.432816;
D = -2.772874;
E = -0.158558;
Cp = A + B*t + C*t^2 + D*t^3 + E/t^2; % J/mol`K
Cp = Cp / MH2 * 1000;
end

function Cp = CpAr(T) % J/kg`K
global MAr
% DOI:10.18434/T42S31
t = T/1000;
A = 20.78600;
B = 2.825911E-7;
C = -1.464191E-7;
D = 1.092131E-8;
E = -3.661371E-8;
Cp = A + B*t + C*t^2 + D*t^3 + E/t^2; % J/mol`K
Cp = Cp / MAr * 1000;
end

function Cp = Cpm(T) % J/kg`K
global ArToH2VolRatio
Cp = (rhoH2(T)*CpH2(T) + ArToH2VolRatio*rhoAr(T)*CpAr(T))/(rhoH2(T)+ArToH2VolRatio*rhoAr(T));
end

function mu = muH2(T) % Pa`s
% DOI:10.1063/1.3160306
mu = -7.118E-12*T^2 + 2.449E-8*T + 2.222E-6;
end

function mu = muAr(T) % Pa`s
% DOI:10.1063/1.556037
mu = -2.608E-11*T^2 + 7.799E-8*T + 1.565E-6;
end

function mu = mum(T)
global MH2 MAr ArToH2VolRatio
% Wilke 1950, DOI:10.1063/1.1747673
phi12 = (1+(muH2(T)/muAr(T)*rhoAr(T)/rhoH2(T))^0.5*(MH2/MAr)^0.25)^2/(4/sqrt(2))/(1+MH2/MAr)^0.5;
phi21 = (1+(muAr(T)/muH2(T)*rhoH2(T)/rhoAr(T))^0.5*(MAr/MH2)^0.25)^2/(4/sqrt(2))/(1+MAr/MH2)^0.5;
mu = muH2(T)/(1+phi12*ArToH2VolRatio) + muAr(T)/(1+phi21/ArToH2VolRatio);
end

function k = kH2(T) % W/m`K
% CRCHCP - Fitted by 2-order polynomial
k = -3.343E-8 * T^2 + 7.092E-4 * T + 0.00224;
end

function k = kAr(T) % W/m`K
% CRCHCP - Fitted by 2-order polynomial
k = -2.839E-8 * T^2 + 6.756E-5 * T - 9E-5;
end

function k = km(T) % W/m`K
global MH2 MAr ArToH2VolRatio TbH2 TbAr
% Lindsay & Bromley 1950, DOI:10.1021/ie50488a017
SH2 = 1.5*TbH2;
SAr = 1.5*TbAr;
S12 = sqrt(SH2*SAr);
A12 = 1/4*(1+(muH2(T)/muAr(T)*(MAr/MH2)^0.75*(1+SH2/T)/(1+SAr/T))^0.5)^2 * (1+S12/T)/(1+SH2/T);
A21 = 1/4*(1+(muAr(T)/muH2(T)*(MH2/MAr)^0.75*(1+SAr/T)/(1+SH2/T))^0.5)^2 * (1+S12/T)/(1+SAr/T);
k = kH2(T)/(1+A12*ArToH2VolRatio) + kAr(T)/(1+A21/ArToH2VolRatio);
end

function nu = num(T)
nu = mum(T)/rhom(T);
end

function alpha = alpham(T)
alpha = km(T)/rhom(T)/Cpm(T);
end

%% = func: Forced Convection Heat Transfer =
function Nu = Nu_fc(Tp) % Nusselt Number, assume laminar flow - Re < 500,000
global l um
vTf = Tf(Tp);
Re = um*l/num(vTf); % Reynolds number
Pr = num(vTf)/alpham(vTf); % Prandtl number
Nu = 0.664*Re^0.5*Pr^0.33;
end

function T = Tf(Tp) % Surface flow temperature = average of Tp and Tm
global Tm
T = (Tp+Tm)/2;
end

%% = func: Natural Convection Heat Transfer =
function Nu = Nu_nc(Tp)
global g Tm Lc flagPlateMaterial
vTf = Tf(Tp);
Pr = num(vTf)/alpham(vTf);
Gr = g*(1/vTf)*(Tp-Tm)*Lc^3/num(Tf)^2;
Nu = 0.54*(Gr*Pr)^0.25; % upper surface
if flagPlateMaterial == 2
    Nu = Nu + 0.27*(Gr*Pr)^0.25; % lower surface 
end
end
