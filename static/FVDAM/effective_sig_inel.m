sig_inel_11=h.*l.*(C11.*D(:,1)+C12.*D(:,2)+C13.*D(:,3) - Rbar(:,1) - (C11.*alpha11+C12.*alpha22+C13.*alpha33)*delta_T);
sig_inel_21=h.*l.*(C12.*D(:,1)+C22.*D(:,2)+C23.*D(:,3) - Rbar(:,2) - (C12.*alpha11+C22.*alpha22+C23.*alpha33)*delta_T);
sig_inel_31=h.*l.*(C13.*D(:,1)+C23.*D(:,2)+C33.*D(:,3) - Rbar(:,3) - (C13.*alpha11+C23.*alpha22+C33.*alpha33)*delta_T);
sig_inel_41=h.*l.*(C44.*D(:,4) - Rbar(:,4));
sig_inel_51=h.*l.*(C55.*D(:,5) - Rbar(:,5));
sig_inel_61=h.*l.*(C66.*D(:,6) - Rbar(:,6));

sig_inel=[sum(sig_inel_11); sum(sig_inel_21); sum(sig_inel_31);....
    sum(sig_inel_41); sum(sig_inel_51); sum(sig_inel_61)];

% T2inv11=T2_inv(1,1);
% T2inv12=T2_inv(1,2);
% T2inv13=T2_inv(1,3);
% T2inv14=T2_inv(1,4);
% T2inv15=T2_inv(1,5);
% T2inv16=T2_inv(1,6);
%
% T2inv21=T2_inv(2,1);
% T2inv22=T2_inv(2,2);
% T2inv23=T2_inv(2,3);
% T2inv24=T2_inv(2,4);
% T2inv25=T2_inv(2,5);
% T2inv26=T2_inv(2,6);
%
% T2inv31=T2_inv(3,1);
% T2inv32=T2_inv(3,2);
% T2inv33=T2_inv(3,3);
% T2inv34=T2_inv(3,4);
% T2inv35=T2_inv(3,5);
% T2inv36=T2_inv(3,6);
%
% T2inv41=T2_inv(4,1);
% T2inv42=T2_inv(4,2);
% T2inv43=T2_inv(4,3);
% T2inv44=T2_inv(4,4);
% T2inv45=T2_inv(4,5);
% T2inv46=T2_inv(4,6);
%
% T2inv51=T2_inv(5,1);
% T2inv52=T2_inv(5,2);
% T2inv53=T2_inv(5,3);
% T2inv54=T2_inv(5,4);
% T2inv55=T2_inv(5,5);
% T2inv56=T2_inv(5,6);
%
% T2inv61=T2_inv(6,1);
% T2inv62=T2_inv(6,2);
% T2inv63=T2_inv(6,3);
% T2inv64=T2_inv(6,4);
% T2inv65=T2_inv(6,5);
% T2inv66=T2_inv(6,6);
%
%
% Dg11=T2inv11*D(:,1)+T2inv12*D(:,2)+T2inv13*D(:,3)+T2inv14*D(:,4)+T2inv15*D(:,5)+T2inv16*D(:,6);
% Dg21=T2inv21*D(:,1)+T2inv22*D(:,2)+T2inv23*D(:,3)+T2inv24*D(:,4)+T2inv25*D(:,5)+T2inv26*D(:,6);
% Dg31=T2inv31*D(:,1)+T2inv32*D(:,2)+T2inv33*D(:,3)+T2inv34*D(:,4)+T2inv35*D(:,5)+T2inv36*D(:,6);
% Dg41=T2inv41*D(:,1)+T2inv42*D(:,2)+T2inv43*D(:,3)+T2inv44*D(:,4)+T2inv45*D(:,5)+T2inv46*D(:,6);
% Dg51=T2inv51*D(:,1)+T2inv52*D(:,2)+T2inv53*D(:,3)+T2inv54*D(:,4)+T2inv55*D(:,5)+T2inv56*D(:,6);
% Dg61=T2inv61*D(:,1)+T2inv62*D(:,2)+T2inv63*D(:,3)+T2inv64*D(:,4)+T2inv65*D(:,5)+T2inv66*D(:,6);
%
% Rbarg11=T2inv11*Rbar(:,1)+T2inv12*Rbar(:,2)+T2inv13*Rbar(:,3)+T2inv14*Rbar(:,4)+T2inv15*Rbar(:,5)+T2inv16*Rbar(:,6);
% Rbarg21=T2inv21*Rbar(:,1)+T2inv22*Rbar(:,2)+T2inv23*Rbar(:,3)+T2inv24*Rbar(:,4)+T2inv25*Rbar(:,5)+T2inv26*Rbar(:,6);
% Rbarg31=T2inv31*Rbar(:,1)+T2inv32*Rbar(:,2)+T2inv33*Rbar(:,3)+T2inv34*Rbar(:,4)+T2inv35*Rbar(:,5)+T2inv36*Rbar(:,6);
% Rbarg41=T2inv41*Rbar(:,1)+T2inv42*Rbar(:,2)+T2inv43*Rbar(:,3)+T2inv44*Rbar(:,4)+T2inv45*Rbar(:,5)+T2inv46*Rbar(:,6);
% Rbarg51=T2inv51*Rbar(:,1)+T2inv52*Rbar(:,2)+T2inv53*Rbar(:,3)+T2inv54*Rbar(:,4)+T2inv55*Rbar(:,5)+T2inv56*Rbar(:,6);
% Rbarg61=T2inv61*Rbar(:,1)+T2inv62*Rbar(:,2)+T2inv63*Rbar(:,3)+T2inv64*Rbar(:,4)+T2inv65*Rbar(:,5)+T2inv66*Rbar(:,6);
%
% alphag11=T2inv11*alpha11+T2inv12*alpha22+T2inv13*alpha33;
% alphag21=T2inv21*alpha11+T2inv22*alpha22+T2inv23*alpha33;
% alphag31=T2inv31*alpha11+T2inv32*alpha22+T2inv33*alpha33;
% alphag41=T2inv41*alpha11+T2inv42*alpha22+T2inv43*alpha33;
% alphag51=T2inv51*alpha11+T2inv52*alpha22+T2inv53*alpha33;
% alphag61=T2inv61*alpha11+T2inv62*alpha22+T2inv63*alpha33;
%
% sig_inel_11=h.*l.*(C11.*Dg11+C12.*Dg21+C13.*Dg31 - Rbarg11 - (C11.*alphag11+C12.*alphag21+C13.*alphag31)*delta_T);
% sig_inel_21=h.*l.*(C12.*Dg11+C22.*Dg21+C23.*Dg31 - Rbarg21 - (C12.*alphag11+C22.*alphag21+C23.*alphag31)*delta_T);
% sig_inel_31=h.*l.*(C13.*Dg11+C23.*Dg21+C33.*Dg31 - Rbarg31 - (C13.*alphag11+C23.*alphag21+C33.*alphag31)*delta_T);
% sig_inel_41=h.*l.*(C44.*Dg41 - Rbarg41 - C44.*alphag41*delta_T);
% sig_inel_51=h.*l.*(C55.*Dg51 - Rbarg51 - C55.*alphag51*delta_T);
% sig_inel_61=h.*l.*(C66.*Dg61 - Rbarg61 - C66.*alphag61*delta_T);
%
% sig_inel=[sum(sig_inel_11); sum(sig_inel_21); sum(sig_inel_31);....
%     sum(sig_inel_41); sum(sig_inel_51); sum(sig_inel_61)];
%


% sig_inel_11=h.*l.*(C11.*D(:,1)+C12.*D(:,2)+C13.*D(:,3)  - (C11.*alpha11+C12.*alpha22+C13.*alpha33)*delta_T);
% sig_inel_21=h.*l.*(C12.*D(:,1)+C22.*D(:,2)+C23.*D(:,3)  - (C12.*alpha11+C22.*alpha22+C23.*alpha33)*delta_T);
% sig_inel_31=h.*l.*(C13.*D(:,1)+C23.*D(:,2)+C33.*D(:,3)  - (C13.*alpha11+C23.*alpha22+C33.*alpha33)*delta_T);
% sig_inel_41=h.*l.*(C44.*D(:,4));
% sig_inel_51=h.*l.*(C55.*D(:,5));
% sig_inel_61=h.*l.*(C66.*D(:,6));
% 
% sig_inel=[sum(sig_inel_11); sum(sig_inel_21); sum(sig_inel_31);....
%     sum(sig_inel_41); sum(sig_inel_51); sum(sig_inel_61)];
% 
% sig_inel = -1/(H*L)*T1_inv*sig_inel;
% 
% Rbar_11=h.*l.*Rbar(:,1);
% Rbar_21=h.*l.*Rbar(:,2);
% Rbar_31=h.*l.*Rbar(:,3);
% Rbar_41=h.*l.*Rbar(:,4);
% Rbar_51=h.*l.*Rbar(:,5);
% Rbar_61=h.*l.*Rbar(:,6);
% 
% Rbar_macro=[sum(Rbar_11); sum(Rbar_21); sum(Rbar_31);....
%     sum(Rbar_41); sum(Rbar_51); sum(Rbar_61)];
% 
% Rbar_global = 1/(H*L)*T2_inv*Rbar_macro;
% 
% sig_inel_bar = sig_inel + Rbar_global;
