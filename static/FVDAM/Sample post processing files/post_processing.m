clear all;

% figure;
% hold on;
% 
% load SiCTi_square_ruc_40_1_offaxis_90_uncracked.mat av*
% plot(100*(av_strain_macro_bar(1,:)), av_stress_macro_bar(1,:),'k-')
% 
% load SiCTi_square_ruc_40_1_offaxis_90_cracked.mat av*
% plot(100*(av_strain_macro_bar(1,:)), av_stress_macro_bar(1,:),'ko','markersize',6)
% 
% load SiCTi_square_ruc_40_1_offaxis_90_cracked_debonded.mat av*
% plot(100*(av_strain_macro_bar(1,:)), av_stress_macro_bar(1,:),'k--','linewidth',1.5)
% 
% load SiCTi_square_ruc_40_1_offaxis_90_uncracked_debonded.mat av*
% plot(100*(av_strain_macro_bar(1,:)), av_stress_macro_bar(1,:),'ko','markerfacecolor',[0 0 0],'markersize',6)
% 
% 
% xlim([0 1.5]);
% xlabel('\bf{\epsilon_{22} (%)}','fontsize',16)
% ylabel('\bf{\sigma_{22} (MPa)}','fontsize',16)
% % legend('Uncracked fiber, perfectly bonded','Cracked fiber, perfectly bonded','Debonding with fiber uncracked',4);
% H=legend('\bf{Uncracked fiber, perfectly bonded}','\bf{Cracked fiber, perfectly bonded}','\bf{Debonding with fiber cracked}','\bf{Debonding with fiber uncracked}',4);
% set(H,'fontsize',14,'box','off');
% % title('Effect of fiber cracking and interface debonding (SiC/Ti Composite)', 'fontsize',12)
% axis('square');
% set(gca,'box','on');
% break;

% figure;
% hold on;
% 
% load SiCTi_square_ruc_40_6_offaxis_90_uncracked.mat av*
% plot(100*(av_strain_macro_bar(6,:)), av_stress_macro_bar(6,:),'k-')
% 
% load SiCTi_square_ruc_40_6_offaxis_90_cracked.mat av*
% plot(100*(av_strain_macro_bar(6,:)), av_stress_macro_bar(6,:),'ko','markersize',6)
% 
% load SiCTi_square_ruc_40_6_offaxis_0_cracked_debonded.mat av*
% plot(100*(av_strain_macro_bar(6,:)), av_stress_macro_bar(6,:),'k--','linewidth',1.5)
% 
% load SiCTi_square_ruc_40_6_offaxis_90_uncracked_debonded.mat av*
% plot(100*(av_strain_macro_bar(6,:)), av_stress_macro_bar(6,:),'ko','markerfacecolor',[0 0 0],'markersize',6)
% 
% 
% xlim([0 1.5]);
% xlabel('\bf{\epsilon_{12} (%)}','fontsize',16)
% ylabel('\bf{\sigma_{12} (MPa)}','fontsize',16)
% % legend('Uncracked fiber, perfectly bonded','Cracked fiber, perfectly bonded','Debonding with fiber cracked','Debonding with fiber uncracked',4);
% % legend('boxoff')
% % title('Effect of fiber cracking and interface debonding (SiC/Ti Composite)', 'fontsize',12)
% axis('square');
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
% plot(100*(av_strain_macro_bar(1,:)), av_stress_macro_bar(1,:),'k^','markerfacecolor',[0 0 0],'markersize',6)
% 
% load SiCTi_square_ruc_40_1_offaxis_45_uncracked.mat av*
% plot(100*(av_strain_macro_bar(1,:)), av_stress_macro_bar(1,:),'ks','markerfacecolor',[0 0 0],'markersize',6)
% 
% load SiCTi_square_ruc_40_1_offaxis_60_uncracked.mat av*
% plot(100*(av_strain_macro_bar(1,:)), av_stress_macro_bar(1,:),'kp','markerfacecolor',[0 0 0],'markersize',6)
% 
% load SiCTi_square_ruc_40_1_offaxis_90_uncracked.mat av*
% plot(100*(av_strain_macro_bar(1,:)), av_stress_macro_bar(1,:),'ko','markerfacecolor',[0 0 0],'markersize',6)
% 
% xlim([0 1.5]);
% xlabel('\bf{\epsilon_{xx} (%)}','fontsize',16)
% ylabel('\bf{\sigma_{xx} (MPa)}','fontsize',16)
% H=legend('\bf{0 \circ}','\bf{10 \circ}','\bf{15 \circ}','\bf{30 \circ}','\bf{45 \circ}','\bf{60 \circ}','\bf{90 \circ}',2);
% set(H,'fontsize',14,'box','off');
% % title('Off-axis response of SiC/Ti composite', 'fontsize',12)
% set(gca,'box','on');
% axis('square');
% % print('temp_plot','-r600','-dtiff');
% break;



% figure;
% hold on;
% 
% load SiCTi_square_ruc_40_1_offaxis_0_uncracked.mat av*
% plot(100*(av_strain_macro_bar(1,:)), av_stress_macro_bar(1,:),'k-')
% 
% load SiCTi_square_ruc_40_1_offaxis_0_cracked.mat av*
% plot(100*(av_strain_macro_bar(1,:)), av_stress_macro_bar(1,:),'ko','markerfacecolor',[0 0 0] ,'markersize',6)
% 
% H=legend('\bf{Uncracked}', '\bf{Cracked}',2);
% set(H,'fontsize',14,'box','off');
% 
% load SiCTi_square_ruc_40_1_offaxis_10_uncracked.mat av*
% plot(100*(av_strain_macro_bar(1,:)), av_stress_macro_bar(1,:),'k-')
% 
% load SiCTi_square_ruc_40_1_offaxis_10_cracked.mat av*
% plot(100*(av_strain_macro_bar(1,:)), av_stress_macro_bar(1,:),'ko','markerfacecolor',[0 0 0] ,'markersize',6)
% 
% load SiCTi_square_ruc_40_1_offaxis_15_uncracked.mat av*
% plot(100*(av_strain_macro_bar(1,:)), av_stress_macro_bar(1,:),'k-')
% 
% load SiCTi_square_ruc_40_1_offaxis_15_cracked.mat av*
% plot(100*(av_strain_macro_bar(1,:)), av_stress_macro_bar(1,:),'ko','markerfacecolor',[0 0 0] ,'markersize',6)
% 
% load SiCTi_square_ruc_40_1_offaxis_30_uncracked.mat av*
% plot(100*(av_strain_macro_bar(1,:)), av_stress_macro_bar(1,:),'k-')
% 
% load SiCTi_square_ruc_40_1_offaxis_30_cracked.mat av*
% plot(100*(av_strain_macro_bar(1,:)), av_stress_macro_bar(1,:),'ko','markerfacecolor',[0 0 0] ,'markersize',6)
% 
% xlim([0 1.5]);
% xlabel('\bf{\epsilon_{xx} (%)}','fontsize',16)
% ylabel('\bf{\sigma_{xx} (MPa)}','fontsize',16)
% axis('square');
% % title('Effect of fiber cracking on the off-axis response of SiC/Ti composite', 'fontsize',12)
% set(gca,'box','on');
% % print('temp_plot','-r600','-dtiff');
% text(1.1,2450,'\bf{0 \circ\rightarrow}','fontsize',14);
% text(1.3,2100,'\bf{\leftarrow{ 10 \circ}}','fontsize',14);
% text(1.2,1650,'\bf{\leftarrow{ 15 \circ}}','fontsize',14);
% text(1.3,950,'\bf{\leftarrow{ 30 \circ}}','fontsize',14);
% break;

figure;
hold on;

load SiCTi_square_ruc_40_1_offaxis_45_uncracked.mat av*
plot(100*(av_strain_macro_bar(1,:)), av_stress_macro_bar(1,:),'k-')

load SiCTi_square_ruc_40_1_offaxis_45_cracked.mat av*
plot(100*(av_strain_macro_bar(1,:)), av_stress_macro_bar(1,:),'k^','markerfacecolor',[0 0 0] ,'markersize',6)

load SiCTi_square_ruc_40_1_offaxis_60_uncracked.mat av*
plot(100*(av_strain_macro_bar(1,:)), av_stress_macro_bar(1,:),'k-.')

load SiCTi_square_ruc_40_1_offaxis_60_cracked.mat av*
plot(100*(av_strain_macro_bar(1,:)), av_stress_macro_bar(1,:),'ks','markerfacecolor',[0 0 0] ,'markersize',6)

load SiCTi_square_ruc_40_1_offaxis_90_uncracked.mat av*
plot(100*(av_strain_macro_bar(1,:)), av_stress_macro_bar(1,:),'k--','linewidth',1.8)

load SiCTi_square_ruc_40_1_offaxis_90_cracked.mat av*
plot(100*(av_strain_macro_bar(1,:)), av_stress_macro_bar(1,:),'ko','markerfacecolor',[0 0 0] ,'markersize',6)

xlim([0 1.5]);
xlabel('\bf{\epsilon_{xx} (%)}','fontsize',16)
ylabel('\bf{\sigma_{xx} (MPa)}','fontsize',16)
H=legend('\bf{45 \circ}','\bf{45 \circ , cracked}','\bf{60 \circ}','\bf{60 \circ , cracked}','\bf{90 \circ}','\bf{90 \circ , cracked}',4);
set(H,'fontsize',14,'box','off');
% title('Effect of fiber cracking on the off-axis response of SiC/Ti composite', 'fontsize',12)
set(gca,'box','on');
axis('square');
% print('temp_plot','-r600','-dtiff');
break;

% load SiCTi_square_ruc_40_3_offaxis_0_uncracked
% figure; plot(av_strain_macro(3,:),av_stress_macro(3,:),'k-','linewidth',1.5); xlabel('\bf{\epsilon_{33}}','fontsize',16); ylabel('\bf{\sigma_{33}}','fontsize',16); 
% break;
% set(figure,'PaperPosition',[0.0 0.0 5 5],'units','inches');
% set(gcf,'Renderer','opengl');
% pcolor(x2,x3,(field_at_incr(25).Sigma22)/abs(av_stress_macro(3,26))), shading interp; 
% % caxis([-500 1800]); 
% caxis([-1 2]); 
% colorbar;
% axis('image');
% axis('off');
% break

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

load SiCTi_square_ruc_40_6_offaxis_90_cracked
% subplot(3,2,1); 
pcolor(x2,x3,field_at_incr(50).Sigma13), shading interp; 
% caxis([0 1500]); 
caxis([-800 800]); 
% colorbar;
axis('image');
axis('off');
% title('\bf{No Cracks, Fully-bonded}');
% print('SiCTi_pressure','-r300','-dtiff');
% close;
% xlabel('\bf{x_2}','fontsize',14)
% ylabel('\bf{x_3}','fontsize',14)
break;

% load SiCTi_square_ruc_40_1_offaxis_90_cracked
% subplot(3,2,3); 
% pcolor(x2,x3,field_at_incr(50).Sigma23), shading interp; 
% % caxis([0 1500]); 
% caxis([-800 800]); 
% caxis([0 0.3]); 
% colorbar;
% axis('image');
% axis('off');
% title('\bf{Radial Fiber Cracks, Fully-bonded}');
% break;

% load SiCTi_square_ruc_40_1_offaxis_90_uncracked_debonded
% subplot(3,2,5); 
pcolor(x2,x3,field_at_incr(50).Sigma23), shading interp; 
% caxis([0 1500]); 
caxis([-800 800]); 
% caxis([0 0.3]); 
colorbar('horiz');
axis('image');
axis('off');
% title('\bf{Debonded, No Cracks}');
