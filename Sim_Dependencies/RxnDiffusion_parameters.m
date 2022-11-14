%%% Written by Agne Tilunaite 2019
%%% Some comments added by Hilary Hunt 2020
function par=RD_parameters(cc0,csr0)

%%% Cytosolic Buffer Parameters
% all [itotal]=uM
% all [kon]=(uMms)^-1
% all [koff]=ms^-1
% Kd=koff/kon

par.ATPtotal=455;   % 455 from Walker, 5000 from Cannell
par.konATP=0.225;   % 0.225 from Walker, 0.1 from Cannell
par.koffATP=45;     % 45 from Walker, 30 from Cannell

par.CaMtotal=24; % 24 from Walker, 36 from Cannell
par.konCaM=0.023;% 0.023 from Walker, 0.1 from Cannell
par.koffCaM=0.238; % 0.238 from Walker, 0.031 from Cannell

par.TnCtotal=70;    % 70 from Walker and Cannell
par.konTnC=0.039;   % 0.039 from Walker, 0.125 from Cannell
par.koffTnC=0.02;   % 0.02 from Walker, 0.25 from Cannell

par.F4total=100;     % 50 from Walker, 100 from Cannell
par.konF4=0.0488;   % 0.1 from Walker, 0.0488 from Cannell
par.koffF4=0.0439;  % 0.11 from Walker, 0.0439 from Cannell

% Luminal Buffer Parameters 
par.CSQtotal=30000; % 30000 from Walker and Cannell
par.konCSQ=0.1;     % 0.1 from Picht and Walker, 0.039 from Cannell
par.koffCSQ=63.8;     % 65 from Picht, 63.8 from Walker, 78 from Cannell
par.KCSQ=800;

par.F5total=25;     % 25 from Hake
par.konF5=250e-6;   % 250e-6 from Hake
par.koffF5=100e-3;  % 100e-3 from Hake

%%% Miscellaneous Buffer Parameters (unused in sim)
par.konSRbuf=0.115;
par.koffSRbuf=0.1;
par.SRbuftotal=0;
par.numserca=1;

par.koffSLbuf=1;
par.konSLbuf=0.115;
par.SLbuftotal=0;

%%% SERCA Parameters (from Walker (cited Tran et al), unused in sim)
par.Kdi=910; %uM
par.Kdsr=2240; %uM
par.Ap=75; %uM
Kij=(cc0/par.Kdi).^2;
Ksr=(csr0/par.Kdsr).^2;
vcycle=(3.24873*(10^12)*(Kij.^2)+Kij.*(9.17846*(10^6)-11478.2.*Ksr)-0.329904.*Ksr)./...
    (0.104217+17.293.*Ksr+Kij.*(1.75583.*(10^6)+7.61673.*(10^6).*Ksr)+(Kij.^2).*(10^11).*(6.08462+4.50544.*Ksr));
par.JSERCA=2*vcycle*par.Ap;

% all diffusivity units in um^2/ms
par.Dcyto=0.22;     % 0.22 from Hake, 0.35 from Cannell
par.Ddyad=0.22;     % 0.22 from Hake, 0.35 from Cannell
par.Dnsr=par.Ddyad/3; % par.Ddyad/3 fitted from Hake et al, 0.06 from Cannell
par.Djsr=0.35;      % 0.35 from Cannell et al
par.DF4Ca=0.042;    % 0.042 from Hake, 0.033 from Cannell
par.DF4=0.042;      % 0.042 from Hake, 0.033 from Cannell
par.DCaM=0.025;     % 0.025 from Hake, 0.045 from Cannell
par.DCaCaM=0.025;   % 0.025 from Hake, 0.045 from Cannell
par.DATP=0.14;      % 0.14 from Hake, 0.15 from Cannell
par.DCaATP=0.14;    % 0.14 from Hake, 0.15 from Cannell
par.DF5=0.008;      % 0.008 from Hake