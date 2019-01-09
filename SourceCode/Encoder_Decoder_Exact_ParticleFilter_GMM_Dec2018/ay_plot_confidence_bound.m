F  = COHB;
G  = COHB_F;
H  = Extra;

ind_g = 12;


subplot(1,2,2)
mf=H.Xf(:,ind_g);
sf=sqrt(H.Wf(:,ind_g));
plot(H.S(15)+exp(mf+2*sf));
plot(H.S(15)+exp(mf-2*sf));
hold off
ay_plot_bound(1,(1:90),(G(:,ind_g))',(H.S(ind_g)+exp(mf-2*sf))',(H.S(ind_g)+exp(mf+2*sf))')
axis tight
hold on
plot(F(:,ind_g),'LineWidth',2)
xlabel('Time')
ylabel('G_c')

subplot(1,2,1)
mf=H.Xf(:,ind_g);
sf=sqrt(H.Wf(:,ind_g));
hold off
ay_plot_bound(1,(1:90),(mf)',(mf-2*sf)',(mf+2*sf)')
axis tight
xlabel('Time')
ylabel('X_{k|K}')




