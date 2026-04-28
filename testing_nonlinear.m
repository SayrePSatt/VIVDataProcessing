iii = 1;
A_y_star = 0.1:0.1:1;
d_sph = 0.08;
m = 2.4564;
m_star = 6.659;
k_eq = 0;
U = 0.1:0.05:0.5;
C_A=0.1;
rho = 998;
delta = 0.39;
f_nw = sqrt(k_eq/m)/(2*pi);
t_span = [0 100];

load lift_fit_4D_obj.mat
lift_fit = lift_fit_4d_obj;
figure
hold on
for ii=1:numel(U)
    for jj=1:numel(A_y_star)
        y_0_nlfreq = [A_y_star(jj)*d_sph(1) 0]; %Displacing from the A_rms value
        [t_nlfreq, disp_nlfreq] = ode45(@(t,y) nonlinear_spring_solver(t,y,m_star,C_A,k_eq(1),lift_fit,d_sph(1),U(ii),rho),t_span,y_0_nlfreq);
        [nl_pks, nl_pks_idx] = findpeaks(disp_nlfreq(:,1));
        f_nl = (numel(nl_pks_idx)-1)/(t_nlfreq(nl_pks_idx(end))-t_nlfreq(nl_pks_idx(1)));
        f_DU(ii,jj) = f_nl*d_sph/U(ii);%*d_sph/U(ii);
        f_w(ii,jj) = sqrt((delta*0.5*rho*U(ii)^2*d_sph)/((m_star+C_A)*rho*4/3*pi*(d_sph/2)^3))/(2*pi)*d_sph/U(ii);
    end
    plot(A_y_star, f_DU(ii,:),"DisplayName",num2str(U(ii)))
end

f_w_line = sqrt(3*delta/(pi*(m_star+C_A)))/(2*pi);
yline(f_w_line)
xlabel('A_rms')
ylabel('fD/U')
% figure
% scatter(f_DU(:),f_w(:))