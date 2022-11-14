%% Ca-dependent parameters of IP3R gating
function [m42_new,h42_new,m24_new,h24_new]=Update_IP3R_vars(Ca_cyto,m42,h42,m24,h24,lambda_h42,ip3r_par,par)

m42_inf=Ca_cyto.^ip3r_par.n42./(ip3r_par.k42.^ip3r_par.n42+Ca_cyto.^ip3r_par.n42);
m42_new=m42+par.dt*ip3r_par.lambda_m42.*(m42_inf-m42);

h42_inf=ip3r_par.kn42.^ip3r_par.nn42./(ip3r_par.kn42^ip3r_par.nn42+Ca_cyto.^ip3r_par.nn42);
h42_new=h42+par.dt*lambda_h42.*(h42_inf-h42); % lambda_h42 is a variable, other lambdas are constants

m24_inf=Ca_cyto.^ip3r_par.n24./(ip3r_par.k24^ip3r_par.n42+Ca_cyto.^ip3r_par.n24);
m24_new=m24+par.dt*ip3r_par.lambda_m24.*(m24_inf-m24);

h24_inf=ip3r_par.kn24.^ip3r_par.nn24./(ip3r_par.kn24.^ip3r_par.nn24+Ca_cyto.^ip3r_par.nn24);
h24_new=h24+par.dt*ip3r_par.lambda_h24.*(h24_inf-h24);

end