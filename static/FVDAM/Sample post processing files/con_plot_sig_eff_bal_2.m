clear all;
% opengl neverselect;

p=figure;
hold on;
colormap(jet);

load('ms_36_SiCTi_square_2_offaxis_0_uncracked');
subplot(1,2,1)
pcolor(x2,x3,field_at_incr(25).Sigma_eff), shading flat;
caxis([400 1000]);  
% axis([min(x2) max(x2) min(x3) max(x3)]);
axis image;
axis 'off'
% G=colorbar; 
% title('\sigma_{12}');

q=figure;
hold on;
colormap(jet);

subplot(1,2,1)
pcolor(x2,x3,field_at_incr(25).Sigma_HS), shading flat;
caxis([0 1200]);  
% axis([min(x2) max(x2) min(x3) max(x3)]);
axis image;
axis 'off'
% H=colorbar; 
% title('\sigma_{12}');

r=figure;
hold on;
colormap(hot);

subplot(1,2,1)
pcolor(x2,x3,field_at_incr(25).Strain_p_eff), shading flat;
caxis([0 0.07]);  
% axis([min(x2) max(x2) min(x3) max(x3)]);
axis image;
axis 'off'
% H=colorbar; 
% title('\sigma_{12}');
clear x* f* l* av*

load('ms_36_SiCTi_circular_2_offaxis_0_uncracked');

figure(p);
subplot(1,2,2)
pcolor(x2,x3,field_at_incr(25).Sigma_eff), shading flat;
caxis([400 1000]);  
% axis([0.25 0.75 0.25 0.75]);
axis image;
axis 'off'
I=colorbar; 
% title('\sigma_{12}');

figure(q);
subplot(1,2,2)
pcolor(x2,x3,field_at_incr(25).Sigma_HS), shading flat;
caxis([0 1200]);  
% axis([0.25 0.75 0.25 0.75]);
axis image;
axis 'off'
J=colorbar; 
% title('\sigma_{12}');

figure(r);
subplot(1,2,2)
pcolor(x2,x3,field_at_incr(25).Strain_p_eff), shading flat;
caxis([0 0.07]);  
% axis([min(x2) max(x2) min(x3) max(x3)]);
axis image;
axis 'off'
H=colorbar; 
% title('\sigma_{12}');
clear x* f* l* av*
