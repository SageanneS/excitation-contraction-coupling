

tic

KConcentration = 4;                                                        % Set extracellular potassium concentration 
Hz = 1;                                                                    % Set stimulation frequency

% Fiber Geometry
endnode = 113;                                                             % Number of Nodes along Sarcolemma
c = 20;                                                                    % Number of Nodes down TT System
r = 20*power(10,-4);                                                       % Fibre Radius
temp = 37;                                                                 % Temperature (Celsius)

% Data to Collect
nodestocollect = ((endnode-1)/2):endnode;

%% File I/O

fout = fopen(['Vm ' num2str(KConcentration) '_' num2str(Hz) '.txt'],'w');

% Read in Fixed ICs
Ids = zeros(endnode,21);
fpin=fopen(['ics ' num2str(KConcentration) '.txt'],'r');
nis=0;
while ~feof(fpin)
    curr = fscanf(fpin,'%f',1);
    if ~isempty(curr)
        nis=nis+1;
        Ids(1,nis)=curr;
    end
end
fclose(fpin);

%% Parameters

current = 3000;                                                            % stimulus current
dt = .001;                                                                 % time step
samprate = 1.0/dt;                                                         % sampling rate
duration = 10;                                                             % simulation duration

% Storing Initialization for Speed
du = zeros(endnode,c); dVm = zeros(1,endnode); gl = zeros(1,c); 
b =  zeros(1,c); Vol = zeros(1,c); Ait = zeros(1,c); u = zeros(endnode,c);
dCa = zeros(endnode,c); dK = zeros(endnode,c); dKs = zeros(1,endnode);
dNas = zeros(1,endnode); dNa = zeros(endnode,c);

% Temperature Coefficients
T = 273 + temp;

q10gi = 1.37;
q10c = 1.02;                                                           
q10k = 2.5;                                                            
q10Na = 2.3;                                                           
q10Na1 = 1.5;                                                          
q10k1 = 1.5;                                                           
q101cl = 1.5;                                                          
q10ir = 1.55;                                                          
q10pump = 1.0;                                                         

temp_coeffgi = power(q10gi, ((temp - 20.0) / 10.0));
temp_coeffc = power(q10c, ((temp - 20.0) / 10.0));
temp_coeffk = power(q10k, ((temp - 20.0) / 10.0));
temp_coeffk1 = power(q10k1, ((temp - 20.0) / 10.0));
temp_coeffNa = power(q10Na, ((temp - 20.0) / 10.0));
temp_coeffNa1 = power(q10Na1, ((temp - 20.0) / 10.0));
temp_coeffir = power(q10ir, ((temp - 20.0) / 10.0));
temp_coeff1cl = power(q101cl, ((temp - 20.0)/ 10.0));
temp_coeffpump = power(q10pump, ((temp - 20.0)/ 10.0));


% General Parameters
Vrp = Ids(1);                                                              % resting potential
V = Vrp*ones(1,endnode);                                                   % initialise membrane voltage at resting potential
Im = zeros(1,endnode);                                                     % initialise stimulus current as an array of zeros for each node
ko = KConcentration*ones(1,endnode);                                       % initialise extracellular potassium concentration 

Cm = 0.58*temp_coeffc;                                                     % membrane capacitance 
F = 96.485;                                                                % Faradays constant
R = 8.31441;                                                               % gas constant
Ri = .125/temp_coeffgi;                                                    % internal resistance
Ra = .150/temp_coeffgi;                                                    % access resistance
dx = 100*power(10,-4);                                                     % internodal distance 
Vsr = power(10,-6);                                                        % volume-surface ratio t-tubular system
VsrS = 4.1*power(10,-6);                                                   % volume-surface ratio sarcolemma
scale = 1000*F*1000*Vsr;
scale_s = 1000*F*1000*VsrS;
p = .003;                                                                  % fraction of muscle volume occupied by t-tubular system
ot = .34;                                                                  % tortuosity factor
Gl = 3.7*p*ot*temp_coeffgi;                                                % luminal conductivity of t-tubular system
Area = 2*pi*r*dx;                                                          % Area of each compartment 
mr = r/((2*Ri*Cm)*(power(dx,2)));                                          % mesh ratio

% Potassium Delayed Rectifier
ki = 154.5*ones(1,endnode);                                                % intracellular potassium concentration                                                        
Ek = ((R*T)/F)*log(ko./ki);                                                % Potassium Nernst potential
gk = temp_coeffk1*21.6;                                                    % max potassium channel conductance
nk = .45;                                                                  % ratio of delayed rectifer potassium channel in t-tubular system
tk = 559;                                                                  % time constant of potassium diffusion between adjacent t-tubular compartments 

% Potassium Delayed Rectifier Rate equations
% Gating Parameters for n (variable of activation behaviour)
a_n = 0.0131*temp_coeffk; b_n = 0.067*temp_coeffk;
Van = 7.0; Vbn = 40.0; Vn = -40; 
an = (a_n*(V(1)-Vn))/(1-exp(-(V(1)-Vn)/Van));                             
bn = b_n*exp(-(V(1)-Vn)/Vbn);
% Gating Parameters for hk (variable of inactivation behaviour
Vhk = -40.0; Ahk = 7.5;
hkinf = (1/(1+exp((V(1)-Vhk)/Ahk)));

% Potassium Inward Rectifier
GK = 3.7*temp_coeffir;
Kk = 950.0;                                                                % dissociation constants for potassium ions
Ks = 1;                                                                    % dissociation constants for blocking cation S+ 
Si = 10.0;                                                                 % intracellular blocking cation concentration
sigma = 0.4;                                                               % fraction of electrical distance through membrane from outside
nir = 1.0;                                                                 % ratio of inward rectifier potassium channel in t-tubular system
ks =  Ks/(power(Si,2)*exp(2*(1-sigma)*(V(1)*F)/(R*T)));                    % first fraction in y equation in Wallinga Paper
Kr = ko(1)*exp(-(sigma*Ek(1)*F)/(R*T));                                    % [Kr] concentration of potassium at binding site
Kr1 = 1 + (power(Kr,2)/Kk);                                                % second fraction in y equation in Wallinga Paper
gir = (GK*power(Kr,2))/(Kk+power(Kr, 2));                                  % surface maximum conductance of inward rectifier channel
y = 1 - power((1 + (ks*Kr1)), -1);                                         % fraction of inward rectifier channels that are open
Gkir = gir.*y;                                                             % surface conductance of inward rectifier channel

% Sodium
gna = temp_coeffNa1*268;                                                   % surface maximum conductance of sodium channel                                                          
nNa = 0.1;                                                                 % ratio of sodium channel density in t-tubular system
Nai = 14.7;                                                                % intracellular sodium concentration  
Nao = 147;                                                                 % extracellular sodium concentration
Ena = ((R*T)/F)*log(Nao./Nai);                                             % Sodium nernst potential

% Sodium Rate equations
% Gating parameters for m (variable of activation behaviour)
a_m = 0.288*temp_coeffNa; b_m = 1.38*temp_coeffNa;
Vm = -46; Vam = 10; Vbm = 18; 
am = a_m*((V(1)-Vm)/(1-exp(-(V(1)-Vm)/Vam)));
bm = b_m*exp(-(V(1)-Vm)/Vbm);
% Gating parameters for h (variable of fast inactivation behaviour
a_h = 0.0081*temp_coeffNa; b_h = 4.38*temp_coeffNa;
Vh = -45; Vah = 14.7; Vbh = 9.0;                                                             
ah = a_h*exp(-(V(1)-Vh)/Vah);
bh = (b_h)/(1+exp(-(V(1)-Vh)/Vbh));
% Gating Parameters for S (variable of slow inactivation behaviour
Vs = -68; As = 7.1; 
Sinf = 1/(1+exp((V(1)-Vs)/As));

% Chloride
nCl = 0.1;
gcl = temp_coeff1cl*6.55;
Va = 70;
Aa = 150;
A = 1/(1+exp((V(1)-Va)/Aa));
Cli = 5.83;
Clo = 128;
Ecl = -((R*T)/F)*log(Clo./Cli);

% NaK Pump
npump = 0.1;                                                               % ratio of channel density for sodium potassium pump in t-tubular system
Jnak = temp_coeffpump*0.207;                                               % maximum pump activity 
Kmk = 1;                                                                   % affinity constant for potassium
KmNa = 13;                                                                 % affinity constant for sodium
I_nak = (F*Jnak)/((1+((Kmk/ko(1))^2))*(1.0+((KmNa/Nai(1))^3)));            % voltage independent part of sodium potassium pump current
o = (exp(Nao(1)/67.3)-1)*(1.0/7.0);                                        % variable in voltage sensitive part of sodium-potassium pump current dependent on extracellular sodium
fnak = (power((1 + (0.12.*exp(-(.1*F*V(1))/(R*T))) + (.04*o*exp(-(F*V(1))/(R*T)))), -1)); % voltage dependent part of the sodium-potassium pump current

%Rate Equations
n = (an/(an+bn))*ones(1,endnode);
m = (am/(am+bm))*ones(1,endnode);
h = (ah/(ah+bh))*ones(1,endnode);
hk = hkinf*ones(1,endnode);
S = Sinf*ones(1,endnode);
fnak = fnak*ones(1,endnode);

%Surface Conductances
Gk=gk*(n(1)^4)*hk(1);
GCl=gcl*(A^4);

%Surface currents
IKDR = Gk.*(V-Ek);
IKIR = Gkir.*(V-Ek);
Gna = gna*(m(1)^3)*h(1)*S(1);
INA = Gna.*(V-Ena);
ICL = GCl.*(V-Ecl);
Inak = I_nak.*fnak;

I = -Im + IKDR + IKIR + INA + ICL + Inak;

%% *******************************************T-system***************************************************

%Tubular Geometry
Vol(1) = pi*dx*((r/20)^2);                                                 % volume of innermost t-system compartment
Ait(1) = (p*Vol(1))/Vsr;                                                   % amount of tubular membrane surface area in innermost compartment of t-system
gl(1)  = (2*pi*(r/20)*dx*Gl)/(r/20);
b(1)   = gl(1)/Ait(1);
for j=2:20
    Vol(j)=pi*dx*((((j*r)/20)^2)-((((j-1)*r)/20)^2));                      % volume of jth t-system compartment
    Ait(j)=(p*Vol(j))/Vsr;                                                 % amount of tubular membrane surface area in jth t-system compartment
    gl(j)=((2*pi*j*(r/20)*dx*Gl)/(r/20));                                  % radial conductivity
    b(j)=gl(j)/Ait(j);
end

%Define Tubular Variables
Ko_t = ko(1)*ones(endnode,c);                                              % extracellular potassium concentration in jth compartment of t-system at ith node (ith segment of muscle fibre)
Ki_t = ki(1)*ones(endnode,c);                                              % intracellular potassium concentration in jth compartment of t-system at ith node (ith segment of muscle fibre)
Ek_t = ((R*T)/F)*log(Ko_t./Ki_t);                                          % Potassium Nernst potential in jth compartment of t-system at ith node (ith segment of muscle fibre)                    

for i = 1:20    
    u(:,i) = Ids(1,22-i)*ones(endnode,1) ; 
end

% Tubular Rate Equations
ant = (a_n.*(u-Vn))./(1-exp(-(u-Vn)./Van));
bnt = b_n.*exp(-(u-Vn)./Vbn);
nt  = ant./(ant+bnt);
amt = a_m.*((u-Vm)./(1-exp(-(u-Vm)./Vam)));
bmt = b_m.*exp(-(u-Vm)./Vbm);
mt  = amt./(amt+bmt);
aht = a_h.*exp(-(u-Vh)./Vah);
bht = b_h./(1+exp(-(u-Vh)./Vbh));
ht  = aht./(aht+bht);
hkinft = 1./(1+exp((u-Vhk)./Ahk));
hkt = hkinft;
Sinft = 1./(1+exp((u-Vs)./As));
St = Sinft;
ft = 1./((1+0.12.*exp(-((.1*F).*u)./(R*T))+(.04*o).*exp(-(F.*u)./(R*T))));

At = 1./(1+exp((u-Va)./Aa));
Cli_t = Cli;
Clo_t = Clo;
Ecl_t = -((R*T)/F).*log(Clo_t./Cli_t);

Nao_t = Nao(1)*ones(endnode,c);
Nai_t = Nai(1)*ones(endnode,c);
Ena_t = ((R*T)/F).*log(Nao_t./Nai_t);

% Tubular Conductances
gkt = nk*gk;
gnat = nNa*gna; 
gclt = nCl*gcl;
Jnakt = npump*Jnak;
GKt = nir*GK;

% Tubular Currents   
Gkt = gkt.*power(nt, 4).*hkt; 
Gnat = gnat.*power(mt, 3).*ht.*St;
GClt = gclt.*power(At,4);
Krt    = Ko_t(1)*exp(-(sigma*Ek_t(1)*F)/(R*T));
Krt1   = 1+((Krt^2)/Kk);
girt   = (GK*(Krt^2))/(Kk+(Krt^2));
kst    = Ks./((Si.^2).*exp(2.*(1-sigma).*(u.*F)./(R*T)));
yt     = 1-((1+(kst.*Krt1)).^-1);
Gkirt  = girt.*yt;
I_nakt = (F*Jnakt)/((1+((Kmk/Ko_t(1))^2))*(1.0+((KmNa/Nai(1))^3)));        % voltage independent part of sodium potassium pump current

IKDRT = Gkt.*(u - Ek_t);
IKIRT = Gkirt.*(u - Ek_t);
INAT = Gnat.*(u - Ena_t);
ICLT = GClt.*(u - Ecl_t);

%% Calcium Dynamics
Katp2 = 0.04;                                                              % uM ATP Dependence of SERCA Pump
Katp3 = 3.5;                                                               % uM ATP Dependence of NaK Pump
K_h =  1000;                                                               % uM Michaelis-Menten constant for ATP hydrolysis
k_HYD = 0.1;                                                               % uM/ms Maximal rate of ATP hydrolysis

% DHPR GATING
kL     = .002;                                                             % /ms Rate Constant for RyR Channel Opening
kLM    = 1000;                                                             % /ms Rate Constant for RyR Channel Closing
fallo  = .2;                                                               % Allosteric Factor         
alpha1 = .2;                                                               % /ms Open RyR Channel Activation Rate Constant
alpha2 = .2;                                                               % /ms Close RyR Channel Activation Rate Constant
KRyR   = 4.5;                                                              % mV RyR Channel Activation Rate Constant
Vbar   = -20;                                                              % mV RyR Channel Activation Rate Constant
i2     = 60.0;                                                             % um^3/ms Fast Twitch RyR Ca Current

% GEOMETRY
Lx       = 1.1;                                                            % um sarcomere length (z-line to m-line)
RR       = 0.5;                                                            % um myofibre radius
V0       = 0.95*(Lx*pi*(RR^2));                                            % um^3 Total myoplasmic volume
VsrC     = 0.05*(Lx*pi*(RR^2));                                            % um^3 SR Volume
V1       = 0.01*V0;                                                        % um^3 TSR Volume     
V2       = 0.99*V0;                                                        % um^3 Myoplasm less TSR less Mitochondria Volume
Vsr1     = 0.01*VsrC;                                                      % um^3 TSR SR Volume 
Vsr2     = 0.99*VsrC;                                                      % um^3 SR less TSR Volume

% CA TRANSPORT AND XB DYNAMICS
Le       = 0.00004;                                                        % uM/ms*um^3 SR Ca leak constant
kCson    = 0.000004;                                                       % /uM*ms Rate of SR Ca binding from Calsequestrin
kCsoff   = 0.005;                                                          % /ms Rate of SR Ca dissociation from Calsequestrin
Cstot    = 31000.0;                                                        % uM Total [Calsequestrin]
vusr     = 2.4375;                                                         % uM/ms* um^3 Rate constant of the SERCA pump
Kcsr     = 1.0;                                                            % uM Dissociation constant of Ca from SERCA
kCatpon  = 0.15;                                                           % /uM*ms Rate of Ca binding to ATP
kCatpoff = 30;                                                             % /ms Rate of Ca dissociation from ATP
kMatpon  = 0.0015;                                                         % /uMms Rate of Mg binding to ATP
kMatpoff = 0.15;                                                           % /ms Rate of Mg dissociation from ATP
tR       = 0.75;                                                           % um^3/ms Intercompartment Ca diffusion parameter
tSR      = tR;                                                             % um^3/ms Intercompartment Ca diffusion parameter
tCa      = tR;
tATP     = 0.375;                                                          % um^3/ms intercompartmental ATP diffusion parameter 
tMg      = 1.5;                                                            % um^3/ms intercompartmental Mg diffusion parameter
Ttot     = 140;                                                            % uM Total [T] Binding Sites 
kTon     = 0.0885;                                                         % /uMms Rate of Ca binding to Troponin
kToff    = 0.115;                                                          % /ms Rate of Ca dissociation from Troponin
k0on     = 0;                                                              % /ms RU activation rate without two Ca bound
k0off    = 0.15;                                                           % /ms RU deactivation raplotte without two Ca bound
kCaon    = 0.15;                                                           % /ms RU activation rate with two Ca bound
kCaoff   = 0.05;                                                           % /ms RU deactivation rate with two Ca bound
f0       = 0.5;                                                            % /ms Rate of XB attachment
fp       = 5;                                                              % /ms Rate of pre-power stroke XB detachment
h0       = 0.08;                                                           % /ms Forward rate of the power stroke
hp       = 0.06;                                                           % /ms Reverse rate of power stroke
g0       = 0.04;                                                           % /ms Rate of post-power stroke XB detachment
bbp      = 0.00000394;                                                     % /ms Rate of myoplasmic phosphate degradation
kp       = 0.00000362;                                                     % um^3/ms Rate of transport of myoplasmic phosphate into the SR
Ap       = 1.0;                                                            % mM^2/ms Rate of phosphate precipitation
Bp       = 0.0001;                                                         % mM/ms Rate of phosphate precipitate solubilization
PP       = 6.0;                                                            % mM^2 phosphate solubility product

% INITIAL CONDITIONS
c0 = ones(endnode,c); 
c1 = zeros(endnode,c); 
c2 = zeros(endnode,c);                                                    
c3 = zeros(endnode,c); 
c4 = zeros(endnode,c); 
o0 = zeros(endnode,c);                                                  
o1 = zeros(endnode,c); 
o2 = zeros(endnode,c); 
o3 = zeros(endnode,c);                                                  
o4 = zeros(endnode,c);                                                                                                 

Ca1 = 0.1*ones(endnode,c);                                                 % TSR Myoplasm [Ca]
Ca2 = 0.1*ones(endnode,c);                                                 % Myoplasm [Ca]
Ca1SR = 1500.0*ones(endnode,c);                                            % TSR [Ca]
Ca2SR = 1500.0*ones(endnode,c);                                            % SR [Ca]
Ca1CS = 16900*ones(endnode,c);                                             % TSR [Ca-Calsequestrin]     
Ca2CS = 16900*ones(endnode,c);                                             % SR [Ca-Calsequestrin]     
Ca1ATP = 0.4*ones(endnode,c);                                              % TSR Myoplasm [Ca-ATP]                           
Ca2ATP = 0.4*ones(endnode,c);                                              % Myoplasm [Ca-ATP]                           
ATP1 = 799.6*ones(endnode,c);                                              % TSR Myoplasm [ATP] 
ATP2 = 799.6*ones(endnode,c);                                              % Myoplasm [ATP]
Mg1ATP = 7200*ones(endnode,c);                                             % TSR Myoplasm [Mg-ATP]              
Mg2ATP = 7200*ones(endnode,c);                                             % Myoplasm [Mg-ATP]              
Mg1 = 1000*ones(endnode,c);                                                % TSR Myoplasm [Mg]
Mg2 = 1000*ones(endnode,c);                                                % Myoplasm [Mg]
Ca2T = 25.0*ones(endnode,c);                                               % Myoplasm [Ca-Troponin]                
Ca2CaT = 3.0*ones(endnode,c);                                              % Myoplasm [Ca-Ca-Troponin]
D0 = 0.8*ones(endnode,c);                                                  % Myoplasm [Detached Activated RU]
D1 = 1.2*ones(endnode,c);                                                  % Myoplasm [Detached Activated RU-Ca]
D2 = 3.0*ones(endnode,c);                                                  % Myoplasm [Detached Activated RU-Ca-Ca]
A1 = 0.23*ones(endnode,c);                                                 % Myoplasm [Attached Pre-Power Stroke XB]
A2 = 0.23*ones(endnode,c);                                                 % Myoplasm [Attached Post-Power Stroke XB]
P = 15*ones(endnode,c);                                                    % Myoplasmic Phosphate                  
PSR = 0.23*ones(endnode,c);                                                % Myoplasmic Phosphate                  
PCSR = 0.23*ones(endnode,c);                                               % Myoplasmic Phosphate                  


%% Ionic Current

% L-type Ca Channel
gca = 3.13;                                                                % maximum conductance of calcium channel
af = 0.1;                                                                  % /ms Rate constant for Ca inactivation 
Kf = 1.0;                                                                  % uM Calcium Sensitivity of Ca inactivation 
fct = 1./(1 + (af./Kf));                                                   % calcium-inactivation gate steady state
Cao_t = 1300*ones(endnode,c);                                              % initial concentration of extracellular calcium
Eca_t = ((R*T)/F)*log(Cao_t./Ca1);                                         % calcium channel nernst potential
ICAT = gca.*(o0+o1+o2+o3+o4).*fct.*(u-Eca_t);  

Inakt = I_nakt.*ft;    
ij = IKDRT + IKIRT + INAT + ICLT + Inakt + ICAT;
It = (V - transpose(u(:,c)))./Ra;

%% Calculations
%loop for Propagation of action potential in both directions from the middle node when a stimulus current is input
skip = (1000)/Hz;
nstep = 1;

for t=1:(duration*samprate)

    if (t >= ceil(0/dt) && t <= ceil(0.4/dt))
       Im((endnode+1)/2) = current;
    end

    V(1)=Vrp;
    V(endnode)=Vrp;
    Vmold = V;
    uold = u;

    %% TRANSVERSE-TUBULAR SYSTEM

    % ION ACCUMULATION
    kotold = Ko_t;
    caotold = Cao_t;
    ext_Ks_old = ko;

    dCa(:,1) = (ICAT(:,1)/scale) - ((caotold(:,1) - caotold(:,2))/tCa);        
    dK(:,1) = (((IKIRT(:,1) + IKDRT(:,1) - (2*Inakt(:,1)))/scale) - ((kotold(:,1) -kotold(:,2))/tk));            
    for j = 2:19
        dCa(:,j) = (ICAT(:,j)/scale) - ((caotold(:,j)-caotold(:,j+1))/tCa) - ((caotold(:,j) - caotold(:,j-1))/tCa);
        dK(:,j) = (((IKIRT(:,j) + IKDRT(:,j) - (2*Inakt(:,j)))/scale) - ((kotold(:,j) - kotold(:,j+1))/tk) - ((kotold(:,j)-kotold(:,j-1))/tk));
    end
    dCa(:,c) = (ICAT(:,c)/scale) - (caotold(:,c) - caotold(:,c-1))/tCa;
    dK(:,c) = (((IKIRT(:,c) +IKDRT(:,c) -(2*Inakt(:,c) ))/scale) - ((kotold(:,c) - ko(:))/tk) - (kotold(:,c) - kotold(:,c-1))/tk);

    dKi = dK.*(15.5/25.0);

    Ek_t = ((R*T)/F).*log(Ko_t./Ki_t);
    Eca_t = ((R*T)/F).*log(Cao_t./Ca1);        

    % RATE EQUATIONS
    ant = (a_n.*(u-Vn))./(1-exp(-(u-Vn)./Van));
    bnt = b_n.*exp(-(u-Vn)./Vbn); 
    amt = a_m.*((u-Vm)./(1-exp(-(u-Vm)./Vam)));
    bmt = b_m.*exp(-(u-Vm)./Vbm);
    aht = a_h.*exp(-(u-Vh)./Vah);
    bht = (b_h)./(1+exp(-(u-Vh)./Vbh));
    At = 1./(1+exp((u-Va)./Aa));
    hkinft = 1./(1+exp((u-Vhk)./Ahk));
    thkt   = exp(-(u+40)./25.75);
    Sinft  = 1./(1+exp((u-Vs)./As));
    tst    = 60./(0.2+(5.65.*power((u+70)./100,2)));
    fctinf = 1./(1 + (Ca1./(Kf)));
    tfct = 1./(af.*(1 + (Ca1/Kf)));
    o_t = (exp(Nao_t./67.3)-1).*(1.0/7.0);                                           
    ft = 1./(1+0.12.*exp(-((.1*F).*u)./(R*T))+(.04*o_t).*exp(-(F.*u)./(R*T)));
    I_nakt = ((F*Jnakt)./((1+((Kmk./Ko_t).^2)).*(1.0+((KmNa./Nai).^3)))).*(ATP2./(Katp3 + ATP2));  

    % GATING VARIABLES
    dnt = ant.*(1-nt)-bnt.*nt; 
    dmt = amt.*(1-mt)-bmt.*mt;
    dht = aht.*(1-ht)-bht.*ht;
    dhkt   = (hkinft-hkt)./(thkt.*(10^3));
    dSt    = (Sinft-St)./(tst.*(10^3));
    dfct = ((fctinf - fct)./(tfct));       

    % CHANNEL CONDUCTANCES
    Gkt   = gkt.*(power(nt,4)).*hkt;
    Krt    = Ko_t.*exp(-(sigma.*Ek_t.*F)./(R*T));
    Krt1   = 1+((Krt.^2)/Kk);
    girt   = (GK.*(Krt.^2))./(Kk+(Krt.^2));
    kst   = Ks./((Si.^2).*exp((2*(1-sigma)).*(u.*F)./(R*T)));
    yt    = 1-(1./(1+(kst.*Krt1)));
    Gkirt = girt.*yt;
    Gnat  = gnat.*(power(mt,3)).*ht.*St;
    GClt  = gclt.*(power(At,4));   

    % IONIC CURRENTS
    IKDRT = Gkt.*(u-Ek_t);
    IKIRT = Gkirt.*(u-Ek_t);
    INAT  = Gnat.*(u-Ena_t);
    ICLT  = GClt.*(u-Ecl_t);
    Inakt = (I_nakt.*ft);
    ICAT = gca.*(o0+o1+o2+o3+o4).*fct.*(u-Eca_t);

    ij = IKDRT+IKIRT+INAT+ICLT+Inakt+ICAT;

    % TT VOLTAGE
    for i = 2:endnode-1
        du(i,1)=(1/(Cm))*(((uold(i,2)-uold(i,1))*b(1))-ij(i,1));
        for j=2:19
            du(i,j)=(1/(Cm))*(((uold(i,j+1)-uold(i,j))*(gl(j)/Ait(j)))-ij(i,j)-(uold(i,j)-uold(i,j-1))*(gl(j-1)/Ait(j)));
        end
        du(i,c)=(1/(Cm))*(((Vmold(i)-uold(i,c))*((2*pi*r*dx)/(Ra*Ait(c))))-ij(i,c)-(uold(i,c)-uold(i,c-1))*(gl(c-1)/Ait(c)));  
    end

    %% SARCOLEMMA

    for i = 2:endnode-1
         dKs(i) = (((IKIR(i)+IKDR(i)-2*Inak(i))/scale_s)-((ext_Ks_old(i)-ext_Ks_old(i+1))/tk + (ext_Ks_old(i)-ext_Ks_old(i-1))/tk + (ext_Ks_old(i)-Ko_t(i,20))/tk));
    end
    dKsi = dKs.*(15.5/25.0);
    Ek = ((R*T)/F)*log(ko./ki);

    % SARCOLEMMA VOLTAGE
    It = (V - transpose(u(:,c)))./Ra;
    for i = 2:endnode - 1
        dVm(i)=(mr*(Vmold(i+1)-2*Vmold(i)+Vmold(i-1)))-((1/Cm)*(I(i)+It(i)));
    end
    Ic = Cm.*(dVm./dt);

    % RATE EQUATIONS
    an = (a_n.*(V-Vn))./(1-exp(-(V-Vn)./Van));
    bn = b_n.*exp(-(V-Vn)./Vbn);
    am = a_m.*((V-Vm)./(1-exp(-(V-Vm)./Vam)));
    bm = b_m.*exp(-(V-Vm)./Vbm);
    ah = a_h.*exp(-(V-Vh)./Vah);
    bh = (b_h)./(1+exp(-(V-Vh)./Vbh));
    A = 1./(1+exp((V-Va)./Aa));   
    hkinf = (1./(1+exp((V-Vhk)./Ahk)));
    thk = exp(-(V+40)./25.75);
    Sinf = (1./(1+exp((V-Vs)./As)));
    ts  = 60./(0.2+(5.65.*(((V+70)./100).^2)));
    o = (exp(Nao./67.3)-1).*(1.0/7.0);                                          
    fnak = 1./(1+0.12.*exp(-((.1*F).*V)./(R*T))+(.04*o).*exp(-(F.*V)./(R*T)));
    I_nak = ((F*Jnak)./((1+((Kmk./ko).^2))*(1.0+((KmNa./Nai).^3))));           

    % GATING VARIABLES
    dn = an.*(1-n)-bn.*n;
    dm = am.*(1-m)-bm.*m;
    dh = ah.*(1-h)-bh.*h;
    dhk = (hkinf-hk)./(thk*(10^3));
    dS  = (Sinf-S)./(ts*(10^3));

    % CHANNEL CONDUCTANCES
    ks   = Ks./((Si.^2).*exp((2*(1-sigma)).*(V.*F)./(R*T)));
    Kr = ko.*exp(-(sigma.*Ek.*F)./(R*T));                                      
    Kr1 = 1 + (power(Kr,2)./Kk);                                               
    gir = (GK.*power(Kr,2))./(Kk+power(Kr, 2));                               
    y    = 1-(1./((1+(ks.*Kr1))));
    Gkir = gir.*y;

    % IONIC CURRENTS
    IKDR = (gk.*power(n,4).*hk).*(V-Ek);
    IKIR = Gkir.*(V-Ek);
    INA  = (gna.*power(m,3).*h.*S).*(V-Ena);
    ICL = gcl.*power(A,4).*(V-Ecl);
    Inak = I_nak.*fnak;
    I = -Im+IKDR+IKIR+INA+ICL+Inak;
    Imem = I+It+Ic;

    %% Calcium Dynamics
    J_SERCA1 = (vusr.*(Ca1./(Ca1+Kcsr))).*(ATP1./(Katp2 + ATP1));  
    J_SERCA2 = (vusr.*(Ca2./(Ca2+Kcsr))).*(ATP2./(Katp2 + ATP2));

    J_HYD1 = (J_SERCA1./2) + k_HYD.*(ATP1./(ATP1 + K_h));
    J_HYD2 = (Inakt./5) + (J_SERCA2./2) + k_HYD.*(ATP2./(ATP2 + K_h)); 

    % DHPR GATING
    kc = 0.5*alpha1.*exp((u-Vbar)/(8*KRyR)); 
    kcm = 0.5*alpha2.*exp(-(u-Vbar)/(8*KRyR));  
    dc0 = -kL.*c0 + kLM.*o0 - 4.*kc.*c0 + kcm.*c1;       
    do0 = kL.*c0 - kLM.*o0 - (4.*kc./fallo).*o0 + fallo.*kcm.*o1;
    dc1 = 4.*kc.*c0 - kcm.*c1 - (kL/fallo).*c1 + kLM*fallo.*o1 - 3.*kc.*c1 + 2.*kcm.*c2;
    do1 = (kL/fallo).*c1 - (kLM*fallo).*o1 + (4.*kc./fallo).*o0 - fallo.*kcm.*o1 - (3.*kc./fallo).*o1 + (2*fallo).*kcm.*o2;
    dc2 = 3.*kc.*c1 - 2.*kcm.*c2 - (kL/(fallo^2)).*c2 + (kLM*(fallo^2)).*o2 - 2.*kc.*c2 + 3.*kcm.*c3;
    do2 = (3.*kc./fallo).*o1 - (2*fallo).*kcm.*o2 + (kL/(fallo^2)).*c2 - (kLM*(fallo^2)).*o2 - (2.*kc./fallo).*o2 + (3*fallo).*kcm.*o3;
    dc3 = 2.*kc.*c2 - 3.*kcm.*c3 - (kL/(fallo^3)).*c3 + (kLM*(fallo^3)).*o3 - kc.*c3 + 4.*kcm.*c4;
    do3 = (kL/(fallo^3)).*c3 - (kLM*(fallo^3)).*o3 + (2.*kc./fallo).*o2 - (3*fallo).*kcm.*o3 - (kc./fallo).*o3 + (4*fallo).*kcm.*o4;
    dc4 = kc.*c3 - 4.*kcm.*c4 - (kL/(fallo^4)).*c4 + (kLM*(fallo^4)).*o4;
    do4 = (kc./fallo).*o3 - (4*fallo).*kcm.*o4 + (kL/(fallo^4)).*c4 - (kLM*(fallo^4)).*o4;

    % CONSERVATION EQUATION
    T0 = Ttot - Ca2T - Ca2CaT - D0 - D1 - D2 - A1 - A2;  

    % CALCIUM RELEASE AND UPTAKE
    J_RYR = (i2.*((o0+o1+o2+o3+o4).*fct).*(Ca1SR - Ca1));

    % CALCIUM CONCENTRATIONS
    dCa1SR = - J_RYR./Vsr1 + J_SERCA1./Vsr1 - Le.*(Ca1SR - Ca1)./Vsr1 - (kCson.*Ca1SR).*(Cstot - Ca1CS) + kCsoff.*Ca1CS - tSR.*(Ca1SR - Ca2SR)./Vsr1;
    dCa1 = ICAT + (Le.*(Ca1SR - Ca1))./V1 + J_RYR./V1  - J_SERCA1./V1 - (kCatpon.*Ca1).*ATP1 + kCatpoff.*Ca1ATP - tR.*(Ca1 - Ca2)./V1;
    dCa2 = Le.*(Ca2SR-Ca2)./V2 - J_SERCA2./V2 + tR.*(Ca1 - Ca2)./V2 - (kCatpon.*Ca2).*ATP2 + kCatpoff.*Ca2ATP - kTon.*Ca2.*T0 + kToff.*Ca2T - kTon.*Ca2.*Ca2T + kToff.*Ca2CaT - kTon.*Ca2.*D0 + kToff.*D1 - kTon.*Ca2.*D1 + kToff.*D2;             

    % BUFFERING EQUATIONS
    dCa1CS =  (kCson.*Ca1SR).*(Cstot - Ca1CS)- kCsoff.*Ca1CS;    
    dCa2CS =  (kCson.*Ca2SR).*(Cstot - Ca2CS) - kCsoff.*Ca2CS;          
    dCa1ATP = (kCatpon.*Ca1).*ATP1 - kCatpoff.*Ca1ATP - tATP.*((Ca1ATP - Ca2ATP))./V1;
    dCa2ATP = (kCatpon.*Ca2).*ATP2 - kCatpoff.*Ca2ATP + tATP.*((Ca1ATP - Ca2ATP))./V2;      
    dMg1ATP = (kMatpon.*Mg1).*ATP1 - kMatpoff.*Mg1ATP - tATP.*((Mg1ATP - Mg2ATP))./V1;
    dMg2ATP = (kMatpon.*Mg2).*ATP2 - kMatpoff.*Mg2ATP + tATP.*((Mg1ATP - Mg2ATP))./V2;        
    dMg1 = -kMatpon.*Mg1.*ATP1 + kMatpoff.*Mg1ATP - tMg.*((Mg1 - Mg2))./V1;   
    dMg2 = -kMatpon.*Mg2.*ATP2 + kMatpoff.*Mg2ATP + tMg.*((Mg1 - Mg2))./V2;              
    dATP1 =  - J_HYD1 - kCatpon.*Ca1.*ATP1 + kCatpoff.*Ca1ATP - kMatpon.*Mg1.*ATP1 + kMatpoff.*Mg1ATP - tATP.*((ATP1- ATP2)./V1);
    dATP2 = - J_HYD2 - kCatpon.*Ca2.*ATP2 + kCatpoff.*Ca2ATP - kMatpon.*Mg2.*ATP2 + kMatpoff.*Mg2ATP + tATP.*((ATP1-ATP2)./V2);          

    % XB DYNAMICS
    dCa2T = (kTon.*Ca2).*T0 - kToff.*Ca2T - (kTon.*Ca2).*Ca2T + kToff.*Ca2CaT - k0on.*Ca2T + k0off.*D1;           
    dCa2CaT = (kTon.*Ca2).*Ca2T - kToff.*Ca2CaT - kCaon.*Ca2CaT + kCaoff.*D2;   
    dD0 = (-kTon.*Ca2).*D0 + kToff.*D1 + k0on.*T0 - k0off.*D0;
    dD1 = (kTon.*Ca2).*D0 - kToff.*D1 + k0on.*Ca2T - k0off.*D1 - (kTon.*Ca2).*D1 + kToff.*D2;
    dD2 = (kTon.*Ca2).*D1 - kToff.*D2 + kCaon.*Ca2CaT -kCaoff.*D2 -f0.*D2 +fp.*A1 +g0.*A2;  
    dA1 = f0.*D2 - fp.*A1 + hp.*A2 -h0.*A1;
    dA2 = -hp.*A2 + h0.*A1 - g0.*A2;  
    dP =  0.001.*J_HYD2 + 0.001.*(h0.*A1 - hp.*A2)-bbp.*P - kp.*((P - PSR)./V2);   

    % METABOLIC FATIGUE
    if PSR.*(0.001).*Ca2SR >= PP
        dPSR = kp.*((P - PSR))./Vsr2 - Ap.*(PSR.*(0.001).*Ca2SR - PP).*(0.001.*PSR.*Ca2SR);
        dPCSR = Ap.*(PSR.*(0.001).*Ca2SR - PP).*(0.001.*PSR.*Ca2SR);
        dCa2SR = J_SERCA2./Vsr2 - Le.*(Ca2SR - Ca2)./Vsr2 + tSR.*(Ca1SR - Ca2SR)./Vsr2 - ((kCson.*Ca2SR).*(Cstot - Ca2CS) - kCsoff.*Ca2CS) - (1000).*(Ap.*(PSR.*(0.001).*Ca2SR - PP).*(0.001).*PSR.*Ca2SR);
    else
        dPSR = kp.*((P - PSR))./Vsr2 + Bp.*PCSR.*(PP -  PSR.*(0.001).*Ca2SR);
        dPCSR = - Bp.*PCSR.*(PP -  PSR.*(0.001).*Ca2SR);
        dCa2SR = J_SERCA2./Vsr2  - Le.*(Ca2SR - Ca2)./Vsr2 + tSR.*(Ca1SR - Ca2SR)./Vsr2 - ((kCson.*Ca2SR).*(Cstot - Ca2CS) - kCsoff.*Ca2CS) - (1000).*(-Bp.*PCSR.*(PP -  PSR.*(0.001).*Ca2SR));
   end

    %% ODE Solving

    % GATING VARIABLES (TRANSVERSE-TUBULES)
    k1_kot = dt.*dK;
    k2_kot = (dt/2).*(dK + k1_kot./2);
    k3_kot = (dt/2).*(dK + k2_kot./2);
    k4_kot = (dt/2).*(dK + k3_kot./2);
    Ko_t = Ko_t + (k1_kot + 2.*k2_kot + 2.*k3_kot + k4_kot)./6;                  

    k1_kit = dt.*dKi;
    k2_kit = (dt/2).*(dKi + k1_kit./2);
    k3_kit = (dt/2).*(dKi + k2_kit./2);
    k4_kit = (dt/2).*(dKi + k3_kit./2);    
    Ki_t = Ki_t - (k1_kit + 2.*k2_kit + 2.*k3_kit + k4_kit)./6;  

    k1_ca = dt.*dCa;
    k2_ca = (dt/2).*(dCa + k1_ca./2);
    k3_ca = (dt/2).*(dCa + k2_ca./2);
    k4_ca = (dt/2).*(dCa + k3_ca./2);
    Cao_t = Cao_t + (k1_ca + 2.*k2_ca + 2.*k3_ca + k4_ca)./6;   

    k1_nt = dt.*dnt;
    k2_nt = (dt/2).*(dnt + k1_nt./2);
    k3_nt = (dt/2).*(dnt + k2_nt./2);
    k4_nt = (dt/2).*(dnt + k3_nt./2);
    nt = nt + (k1_nt + 2.*k2_nt + 2.*k3_nt + k4_nt)./6;

    k1_mt = dt.*dmt;
    k2_mt = (dt/2).*(dmt + k1_mt./2);
    k3_mt = (dt/2).*(dmt + k2_mt./2);
    k4_mt = (dt/2).*(dmt + k3_mt./2);
    mt = mt + (k1_mt + 2.*k2_mt + 2.*k3_mt + k4_mt)./6;

    k1_ht = dt.*dht;
    k2_ht = (dt/2).*(dht + k1_ht./2);
    k3_ht = (dt/2).*(dht + k2_ht./2);
    k4_ht = (dt/2).*(dht + k3_ht./2);
    ht = ht + (k1_ht + 2.*k2_ht + 2.*k3_ht + k4_ht)./6;

    k1_hkt = dt.*dhkt;
    k2_hkt = (dt/2).*(dhkt + k1_hkt./2);
    k3_hkt = (dt/2).*(dhkt + k2_hkt./2);
    k4_hkt = (dt/2).*(dhkt + k3_hkt./2);
    hkt = hkt + (k1_hkt + 2.*k2_hkt + 2.*k3_hkt + k4_hkt)./6;

    k1_St = dt.*dSt;
    k2_St = (dt/2).*(dSt + k1_St./2);
    k3_St = (dt/2).*(dSt + k2_St./2);
    k4_St = (dt/2).*(dSt + k3_St./2);
    St = St + (k1_St + 2.*k2_St + 2.*k3_St + k4_St)./6;

    k1_fct = dt.*dfct;
    k2_fct = (dt/2).*(dfct + k1_fct./2);
    k3_fct = (dt/2).*(dfct + k2_fct./2);
    k4_fct = (dt/2).*(dfct + k3_fct./2);
    fct = fct + (k1_fct + 2.*k2_fct + 2.*k3_fct + k4_fct)./6;

    % GATING VARIABLES (SARCOLEMMA)
    k1_ko = dt.*dKs;
    k2_ko = (dt/2).*(dKs + k1_ko./2);
    k3_ko = (dt/2).*(dKs + k2_ko./2);
    k4_ko = (dt/2).*(dKs + k3_ko./2);
    ko = ko + (k1_ko + 2.*k2_ko + 2.*k3_ko + k4_ko)./6;    

    k1_ki = dt.*dKsi;
    k2_ki = (dt/2).*(dKsi + k1_ki./2);
    k3_ki = (dt/2).*(dKsi + k2_ki./2);
    k4_ki = (dt/2).*(dKsi + k3_ki./2);    
    ki = ki - (k1_ki + 2.*k2_ki + 2.*k3_ki + k4_ki)./6; 

    k1_n = dt.*dn;
    k2_n = (dt/2).*(dn + k1_n./2);
    k3_n = (dt/2).*(dn + k2_n./2);
    k4_n = (dt/2).*(dn + k3_n./2);
    n = n + (k1_n + 2.*k2_n + 2.*k3_n + k4_n)./6; 

    k1_m = dt.*dm;
    k2_m = (dt/2).*(dm + k1_m./2);
    k3_m = (dt/2).*(dm + k2_m./2);
    k4_m = (dt/2).*(dm + k3_m./2);
    m = m + (k1_m + 2.*k2_m + 2.*k3_m + k4_m)./6;

    k1_h = dt.*dh;
    k2_h = (dt/2).*(dh + k1_h./2);
    k3_h = (dt/2).*(dh + k2_h./2);
    k4_h = (dt/2).*(dh + k3_h./2);
    h = h + (k1_h + 2.*k2_h + 2.*k3_h + k4_h)./6;

    k1_hk = dt.*dhk;
    k2_hk = (dt/2).*(dhk + k1_hk./2);
    k3_hk = (dt/2).*(dhk + k2_hk./2);
    k4_hk = (dt/2).*(dhk + k3_hk./2);
    hk = hk + (k1_hk + 2.*k2_hk + 2.*k3_hk + k4_hk)./6;

    k1_S = dt.*dS;
    k2_S = (dt/2).*(dS + k1_S./2);
    k3_S = (dt/2).*(dS + k2_S./2);
    k4_S = (dt/2).*(dS + k3_S./2);
    S = S + (k1_S + 2.*k2_S + 2.*k3_S + k4_S)./6;

    % TRANSMEMBRANE VOLTAGES
    k1_u = dt.*du;
    k2_u = (dt/2).*(du + k1_u./2);
    k3_u = (dt/2).*(du + k2_u./2);
    k4_u = (dt/2).*(du + k3_u./2);
    u = u + (k1_u + 2.*k2_u + 2.*k3_u + k4_u)./6;

    k1_v = dt.*dVm;
    k2_v = (dt/2).*(dVm + k1_v./2);
    k3_v = (dt/2).*(dVm + k2_v./2);
    k4_v = (dt/2).*(dVm + k3_v./2);
    V = V + (k1_v + 2.*k2_v + 2.*k3_v + k4_v)./6;

    % DHPR MODEL
    k1_c0 = dt.*dc0;
    k2_c0 = (dt/2).*(dc0 + k1_c0./2);
    k3_c0 = (dt/2).*(dc0 + k2_c0./2);
    k4_c0 = (dt/2).*(dc0 + k3_c0./2);
    c0 = c0 + (k1_c0 + 2.*k2_c0 + 2.*k3_c0 + k4_c0)./6;

    k1_o0 = dt.*do0;
    k2_o0 = (dt/2).*(do0 + k1_o0./2);
    k3_o0 = (dt/2).*(do0 + k2_o0./2);
    k4_o0 = (dt/2).*(do0 + k3_o0./2);
    o0 = o0 + (k1_o0 + 2.*k2_o0 + 2.*k3_o0 + k4_o0)./6;

    k1_c1 = dt.*dc1;
    k2_c1 = (dt/2).*(dc1 + k1_c1./2);
    k3_c1 = (dt/2).*(dc1 + k2_c1./2);
    k4_c1 = (dt/2).*(dc1 + k3_c1./2);
    c1 = c1 + (k1_c1 + 2.*k2_c1 + 2.*k3_c1 + k4_c1)./6;

    k1_o1 = dt.*do1;
    k2_o1 = (dt/2).*(do1 + k1_o1./2);
    k3_o1 = (dt/2).*(do1 + k2_o1./2);
    k4_o1 = (dt/2).*(do1 + k3_o1./2);
    o1 = o1 + (k1_o1 + 2.*k2_o1 + 2.*k3_o1 + k4_o1)./6;

    k1_c2 = dt.*dc2;
    k2_c2 = (dt/2).*(dc2 + k1_c2./2);
    k3_c2 = (dt/2).*(dc2 + k2_c2./2);
    k4_c2 = (dt/2).*(dc2 + k3_c2./2);
    c2 = c2 + (k1_c2 + 2.*k2_c2 + 2.*k3_c2 + k4_c2)./6;

    k1_o2 = dt.*do2;
    k2_o2 = (dt/2).*(do2 + k1_o2./2);
    k3_o2 = (dt/2).*(do2 + k2_o2./2);
    k4_o2 = (dt/2).*(do2 + k3_o2./2);
    o2 = o2 + (k1_o2 + 2.*k2_o2 + 2.*k3_o2 + k4_o2)./6;

    k1_c3 = dt.*dc3;
    k2_c3 = (dt/2).*(dc3 + k1_c3./2);
    k3_c3 = (dt/2).*(dc3 + k2_c3./2);
    k4_c3 = (dt/2).*(dc3 + k3_c3./2);
    c3 = c3 + (k1_c3 + 2.*k2_c3 + 2.*k3_c3 + k4_c3)./6;

    k1_o3 = dt.*do3;
    k2_o3 = (dt/2).*(do3 + k1_o3./2);
    k3_o3 = (dt/2).*(do3 + k2_o3./2);
    k4_o3 = (dt/2).*(do3 + k3_o3./2);
    o3 = o3 + (k1_o3 + 2.*k2_o3 + 2.*k3_o3 + k4_o3)./6;

    k1_c4 = dt.*dc4;
    k2_c4 = (dt/2).*(dc4 + k1_c4./2);
    k3_c4 = (dt/2).*(dc4 + k2_c4./2);
    k4_c4 = (dt/2).*(dc4 + k3_c4./2);
    c4 = c4 + (k1_c4 + 2.*k2_c4 + 2.*k3_c4 + k4_c4)./6;

    k1_o4 = dt.*do4;
    k2_o4 = (dt/2).*(do4 + k1_o4./2);
    k3_o4 = (dt/2).*(do4 + k2_o4./2);
    k4_o4 = (dt/2).*(do4 + k3_o4./2);
    o4 = o4 + (k1_o4 + 2.*k2_o4 + 2.*k3_o4 + k4_o4)./6;

    % CALCIUM CONCENTRATIONS
    k1_Ca1 = dt.*dCa1;
    k2_Ca1 = (dt/2).*(dCa1 + k1_Ca1/2);
    k3_Ca1 = (dt/2).*(dCa1 + k2_Ca1/2);
    k4_Ca1 = (dt/2).*(dCa1 + k3_Ca1/2);
    Ca1 = Ca1 + (k1_Ca1 + 2*k2_Ca1 + 2*k3_Ca1 + k4_Ca1)/6;

    k1_Ca2 = dt.*dCa2;
    k2_Ca2 = (dt/2).*(dCa2 + k1_Ca2/2);
    k3_Ca2 = (dt/2).*(dCa2 + k2_Ca2/2);
    k4_Ca2 = (dt/2).*(dCa2 + k3_Ca2/2);
    Ca2 = Ca2 + (k1_Ca2 + 2*k2_Ca2 + 2*k3_Ca2 + k4_Ca2)/6;

    k1_Ca1SR = dt.*dCa1SR;
    k2_Ca1SR = (dt/2).*(dCa1SR + k1_Ca1SR/2);
    k3_Ca1SR = (dt/2).*(dCa1SR + k2_Ca1SR/2);
    k4_Ca1SR = (dt/2).*(dCa1SR + k3_Ca1SR/2);
    Ca1SR = Ca1SR + (k1_Ca1SR + 2*k2_Ca1SR + 2*k3_Ca1SR + k4_Ca1SR)/6;

    k1_Ca2SR = dt.*dCa2SR;
    k2_Ca2SR = (dt/2).*(dCa2SR + k1_Ca2SR/2);
    k3_Ca2SR = (dt/2).*(dCa2SR + k2_Ca2SR/2);
    k4_Ca2SR = (dt/2).*(dCa2SR + k3_Ca2SR/2);
    Ca2SR = Ca2SR + (k1_Ca2SR + 2*k2_Ca2SR + 2*k3_Ca2SR + k4_Ca2SR)/6;

    % CALCIUM BUFFERING
    k1_Ca1CS = dt.*dCa1CS;
    k2_Ca1CS = (dt/2).*(dCa1CS + k1_Ca1CS/2);
    k3_Ca1CS = (dt/2).*(dCa1CS + k2_Ca1CS/2);
    k4_Ca1CS = (dt/2).*(dCa1CS + k3_Ca1CS/2);
    Ca1CS = Ca1CS + (k1_Ca1CS + 2*k2_Ca1CS + 2*k3_Ca1CS + k4_Ca1CS)/6;

    k1_Ca2CS = dt.*dCa2CS;
    k2_Ca2CS = (dt/2).*(dCa2CS + k1_Ca2CS/2);
    k3_Ca2CS = (dt/2).*(dCa2CS + k2_Ca2CS/2);
    k4_Ca2CS = (dt/2).*(dCa2CS + k3_Ca2CS/2);
    Ca2CS = Ca2CS + (k1_Ca2CS + 2*k2_Ca2CS + 2*k3_Ca2CS + k4_Ca2CS)/6;

    k1_Ca1ATP = dt.*dCa1ATP;
    k2_Ca1ATP = (dt/2).*(dCa1ATP + k1_Ca1ATP/2);
    k3_Ca1ATP = (dt/2).*(dCa1ATP + k2_Ca1ATP/2);
    k4_Ca1ATP = (dt/2).*(dCa1ATP + k3_Ca1ATP/2);
    Ca1ATP = Ca1ATP + (k1_Ca1ATP + 2*k2_Ca1ATP + 2*k3_Ca1ATP + k4_Ca1ATP)/6;

    k1_Ca2ATP = dt.*dCa2ATP;
    k2_Ca2ATP = (dt/2).*(dCa2ATP + k1_Ca2ATP/2);
    k3_Ca2ATP = (dt/2).*(dCa2ATP + k2_Ca2ATP/2);
    k4_Ca2ATP = (dt/2).*(dCa2ATP + k3_Ca2ATP/2);
    Ca2ATP = Ca2ATP + (k1_Ca2ATP + 2*k2_Ca2ATP + 2*k3_Ca2ATP + k4_Ca2ATP)/6;

    k1_Mg1ATP = dt.*dMg1ATP;
    k2_Mg1ATP = (dt/2).*(dMg1ATP + k1_Mg1ATP/2);
    k3_Mg1ATP = (dt/2).*(dMg1ATP + k2_Mg1ATP/2);
    k4_Mg1ATP = (dt/2).*(dMg1ATP + k3_Mg1ATP/2);
    Mg1ATP = Mg1ATP + (k1_Mg1ATP + 2*k2_Mg1ATP + 2*k3_Mg1ATP + k4_Mg1ATP)/6;

    k1_Mg2ATP = dt.*dMg2ATP;
    k2_Mg2ATP = (dt/2).*(dMg2ATP + k1_Mg2ATP/2);
    k3_Mg2ATP = (dt/2).*(dMg2ATP + k2_Mg2ATP/2);
    k4_Mg2ATP = (dt/2).*(dMg2ATP + k3_Mg2ATP/2);
    Mg2ATP = Mg2ATP + (k1_Mg2ATP + 2*k2_Mg2ATP + 2*k3_Mg2ATP + k4_Mg2ATP)/6;

    k1_Mg1 = dt.*dMg1;
    k2_Mg1 = (dt/2).*(dMg1 + k1_Mg1/2);
    k3_Mg1 = (dt/2).*(dMg1 + k2_Mg1/2);
    k4_Mg1 = (dt/2).*(dMg1 + k3_Mg1/2);
    Mg1 = Mg1 + (k1_Mg1 + 2*k2_Mg1 + 2*k3_Mg1 + k4_Mg1)/6;

    k1_Mg2 = dt.*dMg2;
    k2_Mg2 = (dt/2).*(dMg2 + k1_Mg2/2);
    k3_Mg2 = (dt/2).*(dMg2 + k2_Mg2/2);
    k4_Mg2 = (dt/2).*(dMg2 + k3_Mg2/2);
    Mg2 = Mg2 + (k1_Mg2 + 2*k2_Mg2 + 2*k3_Mg2 + k4_Mg2)/6;

    k1_ATP1 = dt.*dATP1;
    k2_ATP1 = (dt/2).*(dATP1 + k1_ATP1/2);
    k3_ATP1 = (dt/2).*(dATP1 + k2_ATP1/2);
    k4_ATP1 = (dt/2).*(dATP1 + k3_ATP1/2);
    ATP1 = ATP1 + (k1_ATP1 + 2*k2_ATP1 + 2*k3_ATP1 + k4_ATP1)/6;

    k1_ATP2 = dt.*dATP2;
    k2_ATP2 = (dt/2).*(dATP2 + k1_ATP2/2);
    k3_ATP2 = (dt/2).*(dATP2 + k2_ATP2/2);
    k4_ATP2 = (dt/2).*(dATP2 + k3_ATP2/2);
    ATP2 = ATP2 + (k1_ATP2 + 2*k2_ATP2 + 2*k3_ATP2 + k4_ATP2)/6;

    % XB DYNAMICS
    k1_Ca2T = dt.*dCa2T;
    k2_Ca2T = (dt/2).*(dCa2T + k1_Ca2T/2);
    k3_Ca2T = (dt/2).*(dCa2T + k2_Ca2T./2);
    k4_Ca2T = (dt/2).*(dCa2T + k3_Ca2T./2);
    Ca2T = Ca2T + (k1_Ca2T + 2.*k2_Ca2T + 2.*k3_Ca2T + k4_Ca2T)./6;

    k1_Ca2CaT = dt.*dCa2CaT;
    k2_Ca2CaT = (dt/2).*(dCa2CaT + k1_Ca2CaT/2);
    k3_Ca2CaT = (dt/2).*(dCa2CaT + k2_Ca2CaT./2);
    k4_Ca2CaT = (dt/2).*(dCa2CaT + k3_Ca2CaT./2);
    Ca2CaT = Ca2CaT + (k1_Ca2CaT + 2.*k2_Ca2CaT + 2.*k3_Ca2CaT + k4_Ca2CaT)./6;

    k1_D0 = dt.*dD0;
    k2_D0 = (dt/2).*(dD0 + k1_D0./2);
    k3_D0 = (dt/2).*(dD0 + k2_D0./2);
    k4_D0 = (dt/2).*(dD0 + k3_D0./2);
    D0 = D0 + (k1_D0 + 2.*k2_D0 + 2.*k3_D0 + k4_D0)./6;    

    k1_D1 = dt.*dD1;
    k2_D1 = (dt/2).*(dD1 + k1_D1./2);
    k3_D1 = (dt/2).*(dD1 + k2_D1./2);
    k4_D1 = (dt/2).*(dD1 + k3_D1./2);
    D1 = D1 + (k1_D1 + 2.*k2_D1 + 2.*k3_D1 + k4_D1)./6;

    k1_D2 = dt.*dD2;
    k2_D2 = (dt/2).*(dD2 + k1_D2./2);
    k3_D2 = (dt/2).*(dD2 + k2_D2./2);
    k4_D2 = (dt/2).*(dD2 + k3_D2./2);
    D2 = D2 + (k1_D2 + 2.*k2_D2 + 2.*k3_D2 + k4_D2)./6;

    k1_A1 = dt.*dA1;
    k2_A1 = (dt/2).*(dA1 + k1_A1./2);
    k3_A1 = (dt/2).*(dA1 + k2_A1./2);
    k4_A1 = (dt/2).*(dA1 + k3_A1./2);
    A1 = A1 + (k1_A1 + 2.*k2_A1 + 2.*k3_A1 + k4_A1)./6;

    k1_A2 = dt.*dA2;
    k2_A2 = (dt/2).*(dA2 + k1_A2./2);
    k3_A2 = (dt/2).*(dA2 + k2_A2./2);
    k4_A2 = (dt/2).*(dA2 + k3_A2./2);
    A2 = A2 + (k1_A2 + 2.*k2_A2 + 2.*k3_A2 + k4_A2)./6;

    k1_P = dt.*dP;
    k2_P = (dt/2).*(dP + k1_P./2);
    k3_P = (dt/2).*(dP + k2_P./2);
    k4_P = (dt/2).*(dP + k3_P./2);
    P = P + (k1_P + 2.*k2_P + 2.*k3_P + k4_P)./6;    

    k1_PSR = dt.*dPSR;
    k2_PSR = (dt/2).*(dPSR + k1_PSR./2);
    k3_PSR = (dt/2).*(dPSR + k2_PSR./2);
    k4_PSR = (dt/2).*(dPSR + k3_PSR./2);
    PSR = PSR + (k1_PSR + 2.*k2_PSR + 2.*k3_PSR + k4_PSR)./6; 

    k1_PCSR = dt.*dPCSR;
    k2_PCSR = (dt/2).*(dPCSR + k1_PCSR./2);
    k3_PCSR = (dt/2).*(dPCSR + k2_PCSR./2);
    k4_PCSR = (dt/2).*(dPCSR + k3_PCSR./2);
    PCSR = PCSR + (k1_PCSR + 2.*k2_PCSR + 2.*k3_PCSR + k4_PCSR)./6; 

    %% DATA OUTPUT
    
    for i=((endnode-1)/2):endnode
        fprintf(fout,'%f\t',V(i));
    end  
    fprintf(fout,'\n'); 

end
fclose(fout); 

toc    
