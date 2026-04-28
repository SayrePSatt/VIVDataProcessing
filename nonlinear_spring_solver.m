function dydt = nonlinear_spring_solver(t,y,C_A,m_star,k,lift_fit,D,U,rho)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% Pull in the 
    dydt = zeros(2,1);
    total_mass = (m_star+C_A)*4/3*pi*(D/2)^3*rho;
    static_lift = -(rho*U^2*pi*D^2/8)*lift_fit(y(1)/D); %Calculates the nonlinear spring force, 0.5*rho*U^2*D*0.39*y(1);
    dydt(1) = y(2);
    dydt(2) = -k*y(1)/total_mass-static_lift/total_mass;
end