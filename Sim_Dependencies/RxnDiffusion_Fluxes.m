function [Y,J]=RxnDiffusion_Fluxes(Y0,par,par_ip3r,RyR_open,IP3R_open)
% Function calculating main fluxes
% Y0=[Cai, Cajsr, Cansr, CaCaM, CaATP, CaF4, CaTnC, CaCSQ];
%      1     2      3      4      5     6      7      8

% Calculate fluxes for buffers - these are in terms of the forward reaction (taking Ca2+ away)
% [J]=uM/ms
JCaM=(par.koffCaM*(Y0(:,4))-par.konCaM*Y0(:,1).*(par.CaMtotal-Y0(:,4)));    %.*(Y0(:,1)>0.1)
JATP=(par.koffATP*(Y0(:,5))-par.konATP*Y0(:,1).*(par.ATPtotal-Y0(:,5)));
JF4= (par.koffF4*(Y0(:,6))-par.konF4*Y0(:,1).*(par.F4total-Y0(:,6)));
JTnC=(par.koffTnC*(Y0(:,7))-par.konTnC*Y0(:,1).*(par.TnCtotal-Y0(:,7))).*(~par.positions);
JCSQ=(par.koffCSQ*(Y0(:,8))-par.konCSQ*Y0(:,2).*(par.CSQtotal-Y0(:,8))).*par.positions;

% Diffusion matrices (no-flux BC) to be used for:
% Species in dyad
[row_Y0,col_Y0]=size(Y0);
J_diff=zeros(row_Y0,col_Y0);
J_diff(2:end-1,:)=(Y0(3:end,:)+Y0(1:end-2,:)-2*Y0(2:end-1,:))./(par.dx^2);
J_diff(1,:)=(Y0(2,:)-Y0(1,:))./(par.dx^2);           % no flux BC
J_diff(end,:)=(Y0(end-1,:)-Y0(end,:))./(par.dx^2);   % no flux BC

% Species in jSR
J_diff_jSR=zeros(row_Y0,col_Y0);
J_diff_jSR(par.Dyad1(2):par.Dyad1(end-1),:)=(Y0(par.Dyad1(3):par.Dyad1(end),:)+Y0(par.Dyad1(1):par.Dyad1(end-2),:)-2*Y0(par.Dyad1(2):par.Dyad1(end-1),:))./(par.dx^2);
J_diff_jSR(par.Dyad1(1),:)=(Y0(par.Dyad1(2),:)-Y0(par.Dyad1(1),:))./(par.dx^2);
J_diff_jSR(par.Dyad1(end),:)=(Y0(par.Dyad1(end-1),:)-Y0(par.Dyad1(end),:))./(par.dx^2);

J_diff_jSR(par.Dyad2(2):par.Dyad2(end-1),:)=(Y0(par.Dyad2(3):par.Dyad2(end),:)+Y0(par.Dyad2(1):par.Dyad2(end-2),:)-2*Y0(par.Dyad2(2):par.Dyad2(end-1),:))./(par.dx^2);
J_diff_jSR(par.Dyad2(1),:)=(Y0(par.Dyad2(2),:)-Y0(par.Dyad2(1),:))./(par.dx^2);
J_diff_jSR(par.Dyad2(end),:)=(Y0(par.Dyad2(end-1),:)-Y0(par.Dyad2(end),:))./(par.dx^2);

J_diff_jSR(par.Dyad3(2):par.Dyad3(end-1),:)=(Y0(par.Dyad3(3):par.Dyad3(end),:)+Y0(par.Dyad3(1):par.Dyad3(end-2),:)-2*Y0(par.Dyad3(2):par.Dyad3(end-1),:))./(par.dx^2);
J_diff_jSR(par.Dyad3(1),:)=(Y0(par.Dyad3(2),:)-Y0(par.Dyad3(1),:))./(par.dx^2);
J_diff_jSR(par.Dyad3(end),:)=(Y0(par.Dyad3(end-1),:)-Y0(par.Dyad3(end),:))./(par.dx^2);

% Walker's (cited Tran et al) sercas:
Kij=(Y0(:,1)./par.Kdi).^2;
Ksr=(Y0(:,3)./(par.Kdsr)).^2;
vcycle=(3.24873*(10^12)*(Kij.^2)+Kij.*(9.17846*(10^6)-11478.2.*Ksr)-0.329904.*Ksr)./...
    (0.104217+17.293.*Ksr+Kij.*(1.75583.*(10^6)+7.61673.*(10^6).*Ksr)+(Kij.^2).*(10^11).*(6.08462+4.50544.*Ksr));
Jserca=par.scale_sercas*(2*vcycle*par.Ap).*par.sercas;

% Leak
Kij_leak=(0.1*ones(row_Y0,1)./par.Kdi).^2;
Ksr_leak=(Y0(:,3)./(par.Kdsr)).^2;
vcycle_leak=(3.24873*(10^12)*(Kij_leak.^2)+Kij_leak.*(9.17846*(10^6)-11478.2.*Ksr_leak)-0.329904.*Ksr_leak)./...
    (0.104217+17.293.*Ksr_leak+Kij_leak.*(1.75583.*(10^6)+7.61673.*(10^6).*Ksr_leak)+(Kij_leak.^2).*(10^11).*(6.08462+4.50544.*Ksr_leak));
Jleak=par.scale_sercas*(2*vcycle_leak*par.Ap).*par.sercas;

% Release flux from RyR + IP3R
Jrel=(par.gryr*RyR_open+par_ip3r.kipr*IP3R_open).*(Y0(:,2)-Y0(:,1));

% Refill flux from 
Jrefill=(Y0(:,3)-Y0(:,2)).*par.grefill.*par.positions;

% Y0=[Cai, Cajsr, Cansr, CaCaM, CaATP, CaF4, CaTnC, CaCSQ];
%      1     2      3      4      5     6      7      8
Y=[par.Ddyad.*J_diff(:,1)+JF4+JCaM+JATP+JTnC-Jserca+Jleak+Jrel,...
    par.Djsr.*J_diff_jSR(:,2)+JCSQ+Jrefill-Jrel,... 
    par.Dnsr.*J_diff(:,3)+(Jserca-Jrefill-Jleak),...
    par.DCaCaM.*J_diff(:,4)-JCaM,...
    par.DCaATP.*J_diff(:,5)-JATP,...
    par.DF4Ca.*J_diff(:,6)-JF4,...
    -JTnC,...
    -JCSQ];

J=[JCaM,...
    JATP,...
    JF4,...
    JTnC,...
    Jrel,...
    JCSQ,...
    Jserca,...
    Jleak,...
    Jrefill];