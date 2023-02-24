function rho_np1 = lax_friedrichs(rho_n, dt, mesh, L, f)

dx = mesh(2)-mesh(1);
N = length(rho_n);
rho_np1 = NaN(size(rho_n));

v = compute_conv_rep(rho_n,L,mesh);

v = v+f;

rho_n = [rho_n(end); rho_n; rho_n(1)];
v = [v(end); v; v(1)];

fun = rho_n.*v;

for i = 2:N+1
   rho_np1(i-1) = 0.5*(rho_n(i-1)+rho_n(i+1))-dt/(2*dx)*(fun(i+1)-fun(i-1));
end