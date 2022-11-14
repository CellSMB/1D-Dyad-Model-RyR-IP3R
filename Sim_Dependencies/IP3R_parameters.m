%%% Written by Agne Tilunaite 2019
%%% Some comments added by Hilary Hunt 2020
% This function outputs the STEADY STATE parameter values to the Cao et al 2014 IP3R model.

function par_ip3r=IP3R_parameters(c0,ip3c)
% IP3R-2 Transition Rates
par_ip3r.q45=4.14/1000;par_ip3r.q54=3420/1000;
par_ip3r.q12=1140/1000;par_ip3r.q21=95.8/1000;
par_ip3r.q23=4.75/1000;par_ip3r.q32=11.9/1000;
par_ip3r.q26=10100/1000;par_ip3r.q62=3270/1000;

% parameter values
par_ip3r.p=ip3c;    % uM
q24_parameters=[200.3,24.1,46.85,54.9,0.0358,9.5,15.9,774.2,33.0,2.69,1.14,1.14,20.7,1.25,100,0.923,54.5];
q42_parameters=[60,745,1,8.6,0.15,5.5,0.8,19000,86.8,11.6,1.7,37.8,1.2,15.1,1.0,24.5,3.4,2.8];

% NEW PARAMETERS
par_ip3r.v24=(q24_parameters(1)+q24_parameters(2).*par_ip3r.p.^q24_parameters(4)./(par_ip3r.p.^q24_parameters(4)+q24_parameters(3)^q24_parameters(4)))./1000;
par_ip3r.k24=q24_parameters(5);
par_ip3r.n24=q24_parameters(6);
par_ip3r.kn24=q24_parameters(7)+q24_parameters(8)./(par_ip3r.p.^q24_parameters(10)+q24_parameters(9)^q24_parameters(10));   % alters the range at which q24 is low
par_ip3r.nn24=q24_parameters(11)+q24_parameters(12).*par_ip3r.p.^q24_parameters(14)./(par_ip3r.p.^0.25+q24_parameters(13)^q24_parameters(14));
par_ip3r.a24=(q24_parameters(15)./(par_ip3r.p.^q24_parameters(17)+q24_parameters(16)^q24_parameters(17)))./1000;

par_ip3r.v42=(q42_parameters(1)+q42_parameters(2)./(par_ip3r.p.^q42_parameters(4)+q42_parameters(3)^q42_parameters(4)))./1000;
par_ip3r.k42=q42_parameters(5);
par_ip3r.n42=q42_parameters(6);
par_ip3r.kn42=q42_parameters(7)+q42_parameters(8)./(par_ip3r.p.^q42_parameters(10)+q42_parameters(9)^q42_parameters(10));
par_ip3r.nn42=q42_parameters(11)+q42_parameters(12)./(par_ip3r.p.^q42_parameters(14)+q42_parameters(13)^q42_parameters(14));
par_ip3r.a42=(q42_parameters(15)+q42_parameters(16)./(par_ip3r.p.^q42_parameters(18)+q42_parameters(17)^q42_parameters(18)))./1000;

% Scaling Factor if simplified to 2-state model
par_ip3r.scale=(par_ip3r.q62/(par_ip3r.q62+par_ip3r.q26));

% This is the steady state values for used to 
par_ip3r.h24_inf=par_ip3r.kn24^par_ip3r.nn24./(c0^par_ip3r.nn24+par_ip3r.kn24^par_ip3r.nn24);
par_ip3r.lambda_h24=40/1000;
par_ip3r.m24_inf=c0^par_ip3r.n24./(c0^par_ip3r.n24+par_ip3r.k24^par_ip3r.n24);
par_ip3r.lambda_m24=100/1000;

par_ip3r.m42_inf=c0^par_ip3r.n42./(c0^par_ip3r.n42+par_ip3r.k42^par_ip3r.n42);
par_ip3r.lambda_m42=100/1000;
par_ip3r.h42_inf=par_ip3r.kn42^par_ip3r.nn42/(c0^par_ip3r.nn42+par_ip3r.kn42^par_ip3r.nn42);

% lambda_h42 is determined by ah42 or vh42 depending on OPEN or CLOSE
par_ip3r.ah42=0.5/1000;
par_ip3r.vh42=100/1000; % try changing to 100 from 20

% calculates steady state q24, q42 (bc using m_inf and h_inf)
par_ip3r.q24=par_ip3r.a24+par_ip3r.v24.*(1-par_ip3r.m24_inf.*par_ip3r.h24_inf);
par_ip3r.q42=par_ip3r.a42+par_ip3r.v42.*(par_ip3r.m42_inf.*par_ip3r.h42_inf);