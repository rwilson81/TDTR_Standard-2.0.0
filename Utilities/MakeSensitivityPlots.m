for ii=1:length(lambda)
    %-------------Specific heat-----------------
    C_temp=C;
    C_temp(ii)=C(ii)*1.01;
    [deltaR_data_C(:,ii),ratio_data_C(:,ii)]=TDTR_REFL(tdelay,TCR,tau_rep,f,lambda,C_temp,h,eta,r_pump,r_probe,A_pump);
    delta_C(:,ii)=ratio_data_C(:,ii)-ratio_data;
    Num=log(ratio_data_C(:,ii))-log(ratio_data);
    Denom=log(C_temp(ii))-log(C(ii));
    S_C(:,ii)=Num/Denom;
    %-------------Thermal Conductivity (ky)----------
    lambda_temp=lambda;
    lambda_temp(ii)=lambda(ii)*1.01;
    eta_temp = eta.*lambda./lambda_temp;
    [deltaR_data_L(:,ii),ratio_data_L(:,ii)]=TDTR_REFL(tdelay,TCR,tau_rep,f,lambda_temp,C,h,eta_temp,r_pump,r_probe,A_pump);
    delta_L(:,ii)=ratio_data_L(:,ii)-ratio_data;
    Num=log(ratio_data_L(:,ii))-log(ratio_data);
    Denom=log(lambda_temp(ii))-log(lambda(ii));
    S_L(:,ii)=Num/Denom;
    %-------------Layer Thickess---------------
    h_temp=h;
    h_temp(ii)=h(ii)*1.01;
    [deltaR_data_h(:,ii),ratio_data_h(:,ii)]=TDTR_REFL(tdelay,TCR,tau_rep,f,lambda,C,h_temp,eta,r_pump,r_probe,A_pump);
    delta_h(:,ii)=ratio_data_h(:,ii)-ratio_data;
    Num=log(ratio_data_h(:,ii))-log(ratio_data);
    Denom=log(h_temp(ii))-log(h(ii));
    S_h(:,ii)=Num/Denom;
    %--------------------------------------------
    %-------------Anisotropy---------------
 
    eta_temp(ii) = eta(ii)*1.01;
    %eta_temp=eta;
    %eta_temp(ii)=eta_temp(ii)*1.01;
    [deltaR_data_eta(:,ii),ratio_data_eta(:,ii)]=TDTR_REFL(tdelay,TCR,tau_rep,f,lambda,C,h,eta_temp,r_pump,r_probe,A_pump);
    delta_eta(:,ii)=ratio_data_eta(:,ii)-ratio_data;
    Num=log(ratio_data_eta(:,ii))-log(ratio_data);
    Denom=log(eta_temp(ii))-log(eta(ii));
    S_eta(:,ii)=Num/Denom;
    %--------------------------------------------
end
% ftemp=f*1.01;
% [deltaR_data_f,ratio_data_f]=TDTR_REFL(tdelay,TCR,tau_rep,ftemp,lambda,C,h,eta,r_pump,r_probe,A_pump);
% delta_f(:,1)=ratio_data_f-ratio_data;
% Num=log(ratio_data_f)-log(ratio_data);
% Denom=log(ftemp)-log(f);
% S_f=Num/Denom;

r_pumptemp=r_pump*1.01;
r_probetemp=r_probe*1.01
[deltaR_data_r_pump,ratio_data_r_pump]=TDTR_REFL(tdelay,TCR,tau_rep,f,lambda,C,h,eta,r_pumptemp,r_probetemp,A_pump);
delta_r_pump(:,1)=ratio_data_r_pump-ratio_data;
Num=log(ratio_data_r_pump)-log(ratio_data);
Denom=log(r_pumptemp)-log(r_pump);
S_r_pump=Num/Denom;

close all
figure(1)
semilogx(tdelay,[S_C],'*','MarkerSize',8)
hold on
semilogx(tdelay,[S_L],'o','MarkerSize',8)
semilogx(tdelay,[S_h],'x','MarkerSize',8)
semilogx(tdelay,[S_r_pump],'-','LineWidth',2)
semilogx(tdelay,[S_eta],'+','MarkerSize',8)
Cplotlab=strcat('C_',int2str((1:length(lambda))'));
Lplotlab=strcat('kz',int2str((1:length(lambda))'));
tplotlab=strcat('h_',int2str((1:length(lambda))'));
etaplotlab=strcat('kx',int2str((1:length(lambda))'));
%fplotlab='f__';
legend([Cplotlab;Lplotlab;tplotlab;'Rpp';etaplotlab])
set(gca,'FontSize',16)
xlabel('td (ps)','Fontsize',16)
ylabel('Ratio Sensitivity (ps)','FontSize',16)

%dlmwrite('sensitivity.txt',[tdelay,S_C,S_L,S_t,S_r_pump,S_r_probe,S_eta]);


% figure(2)
% semilogx(tdelay,ratio_data,'o')
% axis([100e-12 10e-9 0 max(ratio_data)*1.3])
% title('-Vin/Vout')
% 
% figure(3)
% semilogx(tdelay,real(deltaR_data)/TCR,'b-',tdelay,imag(deltaR_data)/TCR,'g-')
% ampl=sqrt(conj(deltaR_data).*deltaR_data);
% axis([100e-12 10e-9 -max(ampl)*1.3/TCR max(ampl)*1.3/TCR])
% title('Vin')
toc
fprintf('Sensitivities calculated\n')