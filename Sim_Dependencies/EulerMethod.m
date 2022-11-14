function [Yn_plus_1,J]=EulerMethod(Yn,par,ip3r_par,RyR_open,IP3R_open)
% Solves Ca signalling ODEs using 4th order Runge-Kutta method.
% Yn - matrix containing system variables at time tn
% par - structure containing all RD parameters of the model
% RyR/IP3R_open - matrix containing n(open RyR/IP3R) at each element

[K1,J]=RxnDiffusion_Fluxes(Yn,par,ip3r_par,RyR_open,IP3R_open);
Yn_plus_1=Yn+par.dt*K1;

end
