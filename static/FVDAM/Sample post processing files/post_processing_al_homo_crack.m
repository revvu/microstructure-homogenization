% clear all;

% load Al_homo_plastic_ruc_13_offaxis_0_cracked
% plane_strain = field_at_incr;
% load Al_homo_plastic_ruc_3_offaxis_0_cracked
% plane_stress = field_at_incr;

set(figure,'PaperPosition',[0.0 0.0 5 5],'units','inches');
set(gcf,'Renderer','opengl');
colormap('hot');
pcolor(x2,x3,plane_strain(6).Strain_p_eff), shading interp; 
% caxis([0 1500]); 
caxis([0 0.001]); 
% colorbar;
axis('image');
axis('off');
print('al_homo_ruc_cracked_pstrain_p06percent','-r300','-dtiff');
% title('\bf{Debonded, No Cracks}');
close;


break;







xlim([0 1.5]);
xlabel('\bf{\epsilon_{22} (%)}','fontsize',12)
ylabel('\bf{\sigma_{22} (MPa)}','fontsize',12)
legend('Uncracked fiber, perfectly bonded','Cracked fiber, perfectly bonded','Debonding with fiber uncracked',4);
% legend('Uncracked fiber, perfectly bonded','Cracked fiber, perfectly bonded','Debonding with fiber cracked','Debonding with fiber uncracked',4);
legend('boxoff')
title('Effect of fiber cracking and interface debonding (SiC/Ti Composite)', 'fontsize',12)
set(gca,'box','on');
break;

% figure;
% hold on;
% 
% load SiCTi_square_ruc_40_6_offaxis_90_uncracked.mat av*
% plot(100*(av_strain_macro_bar(6,:)), av_stress_macro_bar(6,:),'k-')
% 
% load SiCTi_square_ruc_40_6_offaxis_90_cracked.mat av*
% plot(100*(av_strain_macro_bar(6,:)), av_stress_macro_bar(6,:),'k--')
% 
% % load SiCTi_square_ruc_40_6_offaxis_90_cracked_debonded.mat av*
% % plot(100*(av_strain_macro_bar(6,:)), av_stress_macro_bar(6,:),'r--','linewidth',2.5)
% 
% load SiCTi_square_ruc_40_6_offaxis_90_uncracked_debonded.mat av*
% plot(100*(av_strain_macro_bar(6,:)), av_stress_macro_bar(6,:),'b-')
% 
% 
% xlim([0 1.5]);
% xlabel('\bf{\epsilon_{12} (%)}','fontsize',12)
% ylabel('\bf{\sigma_{12} (MPa)}','fontsize',12)
% % legend('Uncracked fiber, perfectly bonded','Cracked fiber, perfectly bonded','Debonding with fiber cracked','Debonding with fiber uncracked',4);
% % legend('boxoff')
% % title('Effect of fiber cracking and interface debonding (SiC/Ti Composite)', 'fontsize',12)
% set(gca,'box','on');
% break;


% figure;
% hold on;
% 
% load SiCTi_square_ruc_40_1_offaxis_0_uncracked.mat av*
% plot(100*(av_strain_macro_bar(1,:)), av_stress_macro_bar(1,:),'k-')
% 
% load SiCTi_square_ruc_40_1_offaxis_10_uncracked.mat av*
% plot(100*(av_strain_macro_bar(1,:)), av_stress_macro_bar(1,:),'k--')
% 
% load SiCTi_square_ruc_40_1_offaxis_15_uncracked.mat av*
% plot(100*(av_strain_macro_bar(1,:)), av_stress_macro_bar(1,:),'k-.')
% 
% load SiCTi_square_ruc_40_1_offaxis_30_uncracked.mat av*
% plot(100*(av_strain_macro_bar(1,:)), av_stress_macro_bar(1,:),'k^','markerfacecolor',[0 0 0],'markersize',3)
% 
% load SiCTi_square_ruc_40_1_offaxis_45_uncracked.mat av*
% plot(100*(av_strain_macro_bar(1,:)), av_stress_macro_bar(1,:),'ks','markerfacecolor',[0 0 0],'markersize',2)
% 
% load SiCTi_square_ruc_40_1_offaxis_60_uncracked.mat av*
% plot(100*(av_strain_macro_bar(1,:)), av_stress_macro_bar(1,:),'kp','markerfacecolor',[0 0 0],'markersize',3)
% 
% load SiCTi_square_ruc_40_1_offaxis_90_uncracked.mat av*
% plot(100*(av_strain_macro_bar(1,:)), av_stress_macro_bar(1,:),'ko','markerfacecolor',[0 0 0],'markersize',2)
% 
% xlim([0 1.5]);
% xlabel('\bf{\epsilon_{11} (%)}','fontsize',12)
% ylabel('\bf{\sigma_{11} (MPa)}','fontsize',12)
% legend('\bf{0 \circ}','\bf{10 \circ}','\bf{15 \circ}','\bf{30 \circ}','\bf{45 \circ}','\bf{60 \circ}','\bf{90 \circ}',2);
% legend('boxoff')
% title('Off-axis response of SiC/Ti composite', 'fontsize',12)
% set(gca,'box','on');
% % print('temp_plot','-r600','-dtiff');


% figure;
% hold on;
% 
% load SiCTi_square_ruc_40_1_offaxis_0_uncracked.mat av*
% plot(100*(av_strain_macro_bar(1,:)), av_stress_macro_bar(1,:),'k-')
% 
% load SiCTi_square_ruc_40_1_offaxis_0_cracked.mat av*
% plot(100*(av_strain_macro_bar(1,:)), av_stress_macro_bar(1,:),'ko','markerfacecolor',[0 0 0] ,'markersize',2)
% 
% legend('\bf{Uncracked}', '\bf{Cracked}',2);
% legend('boxoff')
% 
% load SiCTi_square_ruc_40_1_offaxis_10_uncracked.mat av*
% plot(100*(av_strain_macro_bar(1,:)), av_stress_macro_bar(1,:),'k-')
% 
% load SiCTi_square_ruc_40_1_offaxis_10_cracked.mat av*
% plot(100*(av_strain_macro_bar(1,:)), av_stress_macro_bar(1,:),'ko','markerfacecolor',[0 0 0] ,'markersize',2)
% 
% load SiCTi_square_ruc_40_1_offaxis_15_uncracked.mat av*
% plot(100*(av_strain_macro_bar(1,:)), av_stress_macro_bar(1,:),'k-')
% 
% load SiCTi_square_ruc_40_1_offaxis_15_cracked.mat av*
% plot(100*(av_strain_macro_bar(1,:)), av_stress_macro_bar(1,:),'ko','markerfacecolor',[0 0 0] ,'markersize',2)
% 
% load SiCTi_square_ruc_40_1_offaxis_30_uncracked.mat av*
% plot(100*(av_strain_macro_bar(1,:)), av_stress_macro_bar(1,:),'k-')
% 
% load SiCTi_square_ruc_40_1_offaxis_30_cracked.mat av*
% plot(100*(av_strain_macro_bar(1,:)), av_stress_macro_bar(1,:),'ko','markerfacecolor',[0 0 0] ,'markersize',2)
% 
% xlim([0 2]);
% xlabel('\bf{\epsilon_{11} (%)}','fontsize',12)
% ylabel('\bf{\sigma_{11} (MPa)}','fontsize',12)
% title('Effect of fiber cracking on the off-axis response of SiC/Ti composite', 'fontsize',12)
% set(gca,'box','on');
% % print('temp_plot','-r600','-dtiff');
% text(1.1,2500,'\bf{0 \circ\rightarrow}','fontsize',14);
% text(1.5,2300,'\bf{\leftarrow{ 10 \circ}}','fontsize',14);
% text(1.2,1600,'\bf{\leftarrow{ 15 \circ}}','fontsize',14);
% text(1.3,900,'\bf{\leftarrow{ 30 \circ}}','fontsize',14);

% figure;
% hold on;
% 
% load SiCTi_square_ruc_40_1_offaxis_45_uncracked.mat av*
% plot(100*(av_strain_macro_bar(1,:)), av_stress_macro_bar(1,:),'k-')
% 
% load SiCTi_square_ruc_40_1_offaxis_45_cracked.mat av*
% plot(100*(av_strain_macro_bar(1,:)), av_stress_macro_bar(1,:),'k^','markerfacecolor',[0 0 0] ,'markersize',3)
% 
% load SiCTi_square_ruc_40_1_offaxis_60_uncracked.mat av*
% plot(100*(av_strain_macro_bar(1,:)), av_stress_macro_bar(1,:),'r-')
% 
% load SiCTi_square_ruc_40_1_offaxis_60_cracked.mat av*
% plot(100*(av_strain_macro_bar(1,:)), av_stress_macro_bar(1,:),'ks','markerfacecolor',[0 0 0] ,'markersize',2)
% 
% load SiCTi_square_ruc_40_1_offaxis_90_uncracked.mat av*
% plot(100*(av_strain_macro_bar(1,:)), av_stress_macro_bar(1,:),'b-')
% 
% load SiCTi_square_ruc_40_1_offaxis_90_cracked.mat av*
% plot(100*(av_strain_macro_bar(1,:)), av_stress_macro_bar(1,:),'ko','markerfacecolor',[0 0 0] ,'markersize',2)
% 
% xlim([0 1.5]);
% xlabel('\bf{\epsilon_{11} (%)}','fontsize',12)
% ylabel('\bf{\sigma_{11} (MPa)}','fontsize',12)
% legend('\bf{45 \circ}','\bf{45 \circ , cracked}','\bf{60 \circ}','\bf{60 \circ , cracked}','\bf{90 \circ}','\bf{90 \circ , cracked}',4);
% legend('boxoff')
% title('Effect of fiber cracking on the off-axis response of SiC/Ti composite', 'fontsize',12)
% set(gca,'box','on');
% % print('temp_plot','-r600','-dtiff');

% load SiCTi_square_ruc_40_3_offaxis_0_uncracked
% set(figure,'PaperPosition',[0.0 0.0 5 5],'units','inches');
% pcolor(x2,x3,field_at_incr(25).Sigma22), shading interp; 
% caxis([-500 1800]); 
% colorbar('horiz');
% axis('image');
% axis('off');
% set(gcf,'Renderer','opengl');
% print('fig5p7_colorbar','-r300','-dtiff');
% close;
% % xlabel('\bf{x_2}','fontsize',14)
% % ylabel('\bf{x_3}','fontsize',14)

% % set(figure,'PaperPosition',[0.0 0.0 8.5 11],'units','inches');
% figure; colormap('hot');
% set(gcf,'Renderer','opengl');

% load SiCTi_square_ruc_40_1_offaxis_90_uncracked
% % subplot(3,2,1); 
% pcolor(x2,x3,field_at_incr(50).Sigma22), shading interp; 
% caxis([0 1500]); 
% % colorbar;
% axis('image');
% axis('off');
% % title('\bf{No Cracks, Fully-bonded}');
% % print('SiCTi_pressure','-r300','-dtiff');
% % close;
% % xlabel('\bf{x_2}','fontsize',14)
% % ylabel('\bf{x_3}','fontsize',14)


% load SiCTi_square_ruc_40_1_offaxis_90_cracked
% % subplot(3,2,3); 
% pcolor(x2,x3,field_at_incr(50).Strain_p_eff), shading interp; 
% % caxis([0 1500]); 
% caxis([0 0.3]); 
% colorbar;
% axis('image');
% axis('off');
% % title('\bf{Radial Fiber Cracks, Fully-bonded}');
% break;
% 
% load SiCTi_square_ruc_40_1_offaxis_90_uncracked_debonded
% % subplot(3,2,5); 
% pcolor(x2,x3,field_at_incr(50).Strain_p_eff), shading interp; 
% % caxis([0 1500]); 
% caxis([0 0.3]); 
% colorbar;
% axis('image');
% axis('off');
% % title('\bf{Debonded, No Cracks}');


% set(figure,'PaperPosition',[0.0 0.0 8.5 11],'units','inches');
figure; colormap('jet');
set(gcf,'Renderer','opengl');

% load SiCTi_square_ruc_40_1_offaxis_90_uncracked
% % subplot(3,2,1); 
% pcolor(x2,x3,field_at_incr(50).Sigma22), shading interp; 
% caxis([0 1500]); 
% % colorbar;
% axis('image');
% axis('off');
% % title('\bf{No Cracks, Fully-bonded}');
% % print('SiCTi_pressure','-r300','-dtiff');
% % close;
% % xlabel('\bf{x_2}','fontsize',14)
% % ylabel('\bf{x_3}','fontsize',14)


% load SiCTi_square_ruc_40_6_offaxis_90_uncracked
% subplot(3,2,3); 
pcolor(x2,x3,field_at_incr(50).Sigma22), shading interp; 
% caxis([0 1500]); 
% caxis([0 0.3]); 
colorbar;
axis('image');
axis('off');
% title('\bf{Radial Fiber Cracks, Fully-bonded}');
break;

% load SiCTi_square_ruc_40_1_offaxis_90_uncracked_debonded
% subplot(3,2,5); 
