iii = 1;
A_y_star(iii) = 1;
d_sph = 0.08;
m = 2.5;
k_eq = 5.8;
U = 0.33;
rho = 998;
f_nw = sqrt(k_eq/m)/(2*pi);

load lift_fit_4D_obj.mat

t_span = [0 10/f_nw(1)];
y_0_nlfreq = [A_y_star(iii)*d_sph(1) 0]; %Displacing from the A_rms value
lift_fit = lift_fit_4d_obj;
[t_nlfreq, disp_nlfreq] = ode45(@(t,y) nonlinear_spring_solver(t,y,m(1),k_eq(1),lift_fit,d_sph(1),U,rho),t_span,y_0_nlfreq);

[nl_pks, nl_pks_idx] = findpeaks(disp_nlfreq(:,1));
f_nl = (numel(nl_pks_idx)-1)/(t_nlfreq(nl_pks_idx(end))-t_nlfreq(nl_pks_idx(1)));
f_star = f_nw/f_nl
plot(t_nlfreq, disp_nlfreq(:,1))