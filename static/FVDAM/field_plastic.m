% Generates the actual displacement and stress fields
function [Sigma11, Sigma22, Sigma33, Sigma23, Sigma13, Sigma12, Sigma_HS, Sigma_eff, Strain_p_eff, u1, u2, u3] = field_plastic(m, n, h, l, G1_1, G2_1, G3_1, G1_2, G2_2, G3_2, C11, C22, C33, C12, C13, C23, C44, C55, C66, C22bar, C33bar, C55bar, C66bar, alpha11, alpha22, alpha33, Avu, Avunext, Ave11, Ave22, Ave33, Ave23, Ave13, Ave12, ul, ur, ub, ut, pl_eps, delta_T, array_x)
% function [Sigma_HS, Sigma_eff, Strain_p_eff] = field_plastic(m, n, h, l, G1_1, G2_1, G3_1, G1_2, G2_2, G3_2, C11, C22, C33, C12, C13, C23, C44, C55, C66, C22bar, C33bar, C55bar, C66bar, alpha11, alpha22, alpha33, Avu, Avunext, Ave11, Ave22, Ave33, Ave23, Ave13, Ave12, ul, ur, ub, ut, pl_eps, delta_T, array_x)
% function [Sigma11, Sigma22, Sigma33, Sigma23, Sigma13, Sigma12, Sigma_HS, Sigma_eff, Strain_p_eff] = field_plastic(m, n, h, l, G1_1, G2_1, G3_1, G1_2, G2_2, G3_2, C11, C22, C33, C12, C13, C23, C44, C55, C66, C22bar, C33bar, C55bar, C66bar, alpha11, alpha22, alpha33, Avu, Avunext, Ave11, Ave22, Ave33, Ave23, Ave13, Ave12, ul, ur, ub, ut, pl_eps, delta_T, array_x)
% function [Sigma22, Sigma33, Sigma23, u2, u3] = field_plastic(m, n, h, l, G1_1, G2_1, G3_1, G1_2, G2_2, G3_2, C11, C22, C33, C12, C13, C23, C44, C55, C66, C22bar, C33bar, C55bar, C66bar, alpha11, alpha22, alpha33, Avu, Avunext, Ave11, Ave22, Ave33, Ave23, Ave13, Ave12, ul, ur, ub, ut, pl_eps, delta_T, array_x)

%%_________________________________________________________________________
NS=m*n;
u221=zeros(NS,1);
u222=zeros(NS,1);
u231=zeros(NS,1);
u232=zeros(NS,1);
u321=zeros(NS,1);
u322=zeros(NS,1);
u331=zeros(NS,1);
u332=zeros(NS,1);
u211=zeros(NS,1);
u212=zeros(NS,1);
u311=zeros(NS,1);
u312=zeros(NS,1);

p=ur(NS)-n;
q=ut(NS)-m;

Avu_22=Avu(1:p);
Avu_23=Avu(p+1:2*p);
Avu_32=Avu(2*p+1:2*p+q);
Avu_33=Avu(2*p+q+1:2*p+2*q);
clear Avu

for i=0:n-1
    t=i*m;
    j=t+1:t+m-1;
    u221(j)=Avu_22(ul(j)-i);
    u222(j)=Avu_22(ur(j)-i);
    u231(j)=Avu_23(ul(j)-i);
    u232(j)=Avu_23(ur(j)-i);
    j=t+m;
    u221(j)=Avu_22(ul(j)-i);
    u222(j)=Avu_22(ul(t+1)-i);
    u231(j)=Avu_23(ul(j)-i);
    u232(j)=Avu_23(ul(t+1)-i);
end

for i=0:n-2
    j=i*m+1:i*m+m;
    k=0:m-1;
    u321(j)=Avu_32(ub(j)-k);
    u322(j)=Avu_32(ut(j)-k);
    u331(j)=Avu_33(ub(j)-k);
    u332(j)=Avu_33(ut(j)-k);
end
j=(n-1)*m+1:m*n;
k=0:m-1;
u321(j)=Avu_32(ub(j)-k);
u322(j)=Avu_32(ub(k+1)-k);
u331(j)=Avu_33(ub(j)-k);
u332(j)=Avu_33(ub(k+1)-k);

Avu_21=Avunext(1:p);
Avu_31=Avunext(p+1:p+q);
clear Avunext

for i=0:n-1
    t=i*m;
    j=t+1:t+m-1;
    u211(j)=Avu_21(ul(j)-i);
    u212(j)=Avu_21(ur(j)-i);
    j=t+m;
    u211(j)=Avu_21(ul(j)-i);
    u212(j)=Avu_21(ul(t+1)-i);
end

for i=0:n-2
    j=i*m+1:i*m+m;
    k=0:m-1;
    u311(j)=Avu_31(ub(j)-k);
    u312(j)=Avu_31(ut(j)-k);
end
j=(n-1)*m+1:m*n;
k=0:m-1;
u311(j)=Avu_31(ub(j)-k);
u312(j)=Avu_31(ub(k+1)-k);

% q=ur(m*n)-n;
% for i=0:n-1
%     t=i*m;
%     j=t+1:t+m-1;
%     u221(j)=Avu(ul(j)-i);
%     u222(j)=Avu(ur(j)-i);
%     u231(j)=Avu(q+ul(j)-i);
%     u232(j)=Avu(q+ur(j)-i);
%     j=t+m;
%     u221(j)=Avu(ul(j)-i);
%     u222(j)=Avu(ul(t+1)-i);
%     u231(j)=Avu(q+ul(j)-i);
%     u232(j)=Avu(q+ul(t+1)-i);
% end
% 
% p=2*(ur(m*n)-n);
% q=ut(m*n)-m;
% for i=0:n-2
%     j=i*m+1:i*m+m;
%     k=0:m-1;
%     u321(j)=Avu(p+ub(j)-k);
%     u322(j)=Avu(p+ut(j)-k);
%     u331(j)=Avu(p+q+ub(j)-k);
%     u332(j)=Avu(p+q+ut(j)-k);
% end
% j=(n-1)*m+1:m*n;
% k=0:m-1;
% u321(j)=Avu(p+ub(j)-k);
% u322(j)=Avu(p+ub(k+1)-k);
% u331(j)=Avu(p+q+ub(j)-k);
% u332(j)=Avu(p+q+ub(k+1)-k);
% 
% 
% for i=0:n-1
%     t=i*m;
%     j=t+1:t+m-1;
%     u211(j)=Avunext(ul(j)-i);
%     u212(j)=Avunext(ur(j)-i);
%     j=t+m;
%     u211(j)=Avunext(ul(j)-i);
%     u212(j)=Avunext(ul(t+1)-i);
% end
% 
% q=ur(m*n)-n;
% for i=0:n-2
%     j=i*m+1:i*m+m;
%     k=0:m-1;
%     u311(j)=Avunext(q+ub(j)-k);
%     u312(j)=Avunext(q+ut(j)-k);
% end
% j=(n-1)*m+1:m*n;
% k=0:m-1;
% u311(j)=Avunext(q+ub(j)-k);
% u312(j)=Avunext(q+ub(k+1)-k);
%%_________________________________________________________________________


%%_________________________________________________________________________
% W1000=zeros(NS,1);
% W2000=zeros(NS,1);
% W3000=zeros(NS,1);
% W1010=zeros(NS,1);
% W2010=zeros(NS,1);
% W3010=zeros(NS,1);
% W1001=zeros(NS,1);
% W2001=zeros(NS,1);
% W3001=zeros(NS,1);
% W1020=zeros(NS,1);
% W2020=zeros(NS,1);
% W3020=zeros(NS,1);
% W1002=zeros(NS,1);
% W2002=zeros(NS,1);
% W3002=zeros(NS,1);
% for i=1:NS
%     W1000(i)=C66(i)*(u212(i)+u211(i))/(2*C66bar(i))+C55(i)*(u312(i)+u311(i))/(2*C55bar(i)) - h(i)*(G3_1(i)+h(i)/l(i)*G3_2(i))/(6*C66bar(i));
%     W2000(i)=C22(i)*(u222(i)+u221(i))/(2*C22bar(i))+h(i)^2*C44(i)*(u322(i)+u321(i))/(2*l(i)^2*C22bar(i)) - h(i)*(G1_1(i)+h(i)/l(i)*G1_2(i))/(6*C22bar(i));
%     W3000(i)=C33(i)*(u332(i)+u331(i))/(2*C33bar(i))+l(i)^2*C44(i)*(u232(i)+u231(i))/(2*h(i)^2*C33bar(i)) - l(i)*(l(i)/h(i)*G2_1(i)+G2_2(i))/(6*C33bar(i));
% 
%     W2010(i)=(u222(i)-u221(i))/h(i);
%     W3010(i)=(u232(i)-u231(i))/h(i);
%     W2001(i)=(u322(i)-u321(i))/l(i);
%     W3001(i)=(u332(i)-u331(i))/l(i);
%     W1010(i)=(u212(i)-u211(i))/h(i);
%     W1001(i)=(u312(i)-u311(i))/l(i);
% 
%     W2020(i)=(u222(i)+u221(i))*2/h(i)^2-4*W2000(i)/h(i)^2;
%     W3020(i)=(u232(i)+u231(i))*2/h(i)^2-4*W3000(i)/h(i)^2;
%     W2002(i)=(u322(i)+u321(i))*2/l(i)^2-4*W2000(i)/l(i)^2;
%     W3002(i)=(u332(i)+u331(i))*2/l(i)^2-4*W3000(i)/l(i)^2;
%     W1020(i)=(u212(i)+u211(i))*2/h(i)^2-4*W1000(i)/h(i)^2;
%     W1002(i)=(u312(i)+u311(i))*2/l(i)^2-4*W1000(i)/l(i)^2;
% end

W1000=C66.*(u212+u211)./(2*C66bar)+C55.*(u312+u311)./(2*C55bar) - h.*(G3_1+h./l.*G3_2)./(6*C66bar);
W2000=C22.*(u222+u221)./(2*C22bar)+h.^2.*C44.*(u322+u321)./(2*l.^2.*C22bar) - h.*(G1_1+h./l.*G1_2)./(6*C22bar);
W3000=C33.*(u332+u331)./(2*C33bar)+l.^2.*C44.*(u232+u231)./(2*h.^2.*C33bar) - l.*(l./h.*G2_1+G2_2)./(6*C33bar);

W2010=(u222-u221)./h;
W3010=(u232-u231)./h;
W2001=(u322-u321)./l;
W3001=(u332-u331)./l;
W1010=(u212-u211)./h;
W1001=(u312-u311)./l;

W2020=2*(u222+u221)./h.^2-4*W2000./h.^2;
W3020=2*(u232+u231)./h.^2-4*W3000./h.^2;
W2002=2*(u322+u321)./l.^2-4*W2000./l.^2;
W3002=2*(u332+u331)./l.^2-4*W3000./l.^2;
W1020=2*(u212+u211)./h.^2-4*W1000./h.^2;
W1002=2*(u312+u311)./l.^2-4*W1000./l.^2;
%%_________________________________________________________________________


t=length(array_x);
[strain_inel_11, strain_inel_22, strain_inel_33, strain_inel_23, strain_inel_13, strain_inel_12, strain_inel_eff] = struct2double3d(m, n, pl_eps, t);


disu1=zeros(t,t,NS);
disu2=zeros(t,t,NS);
disu3=zeros(t,t,NS);
Sig11=zeros(t,t,NS);
Sig22=zeros(t,t,NS);
Sig33=zeros(t,t,NS);
Sig23=zeros(t,t,NS);
Sig13=zeros(t,t,NS);
Sig12=zeros(t,t,NS);
u1=[];
u2=[];
u3=[];
Sigma11=[];
Sigma22=[];
Sigma33=[];
Sigma23=[];
Sigma13=[];
Sigma12=[];
Strain_p_eff=[];
for i=1:n
    u1temp=[];
    u2temp=[];
    u3temp=[];
    Sigm11=[];
    Sigm22=[];
    Sigm33=[];
    Sigm23=[];
    Sigm13=[];
    Sigm12=[];
    Str_p_eff=[];
    for j=1:m

        k=(i-1)*m+j;
        x2=array_x*h(k)/2;
        x3=array_x*l(k)/2;
        [X,Y]=meshgrid(x2,x3);

        disu1(:,:,k)=W1000(k)+W1010(k)*X+W1001(k)*Y+W1020(k)*(3*X.^2-h(k)^2/4)/2+W1002(k)*(3*Y.^2-l(k)^2/4)/2;
        disu2(:,:,k)=W2000(k)+W2010(k)*X+W2001(k)*Y+W2020(k)*(3*X.^2-h(k)^2/4)/2+W2002(k)*(3*Y.^2-l(k)^2/4)/2;
        disu3(:,:,k)=W3000(k)+W3010(k)*X+W3001(k)*Y+W3020(k)*(3*X.^2-h(k)^2/4)/2+W3002(k)*(3*Y.^2-l(k)^2/4)/2;
        Sig11(:,:,k)=C11(k)*Ave11+C12(k)*(Ave22+W2010(k)+3*W2020(k)*X)+C13(k)*(Ave33+W3001(k)+3*W3002(k)*Y) - 2*C44(k)*strain_inel_11(:,:,k) - delta_T*(C11(k)*alpha11(k)+C12(k)*alpha22(k)+C13(k)*alpha33(k));
        Sig22(:,:,k)=C12(k)*Ave11+C22(k)*(Ave22+W2010(k)+3*W2020(k)*X)+C23(k)*(Ave33+W3001(k)+3*W3002(k)*Y) - 2*C44(k)*strain_inel_22(:,:,k) - delta_T*(C12(k)*alpha11(k)+C22(k)*alpha22(k)+C23(k)*alpha33(k));
        Sig33(:,:,k)=C13(k)*Ave11+C23(k)*(Ave22+W2010(k)+3*W2020(k)*X)+C33(k)*(Ave33+W3001(k)+3*W3002(k)*Y) - 2*C44(k)*strain_inel_33(:,:,k) - delta_T*(C13(k)*alpha11(k)+C23(k)*alpha22(k)+C33(k)*alpha33(k));
        Sig23(:,:,k)=C44(k)*(2*Ave23+W2001(k)+W3010(k)+3*W3020(k)*X+3*W2002(k)*Y) - 2*C44(k)*strain_inel_23(:,:,k);
        Sig13(:,:,k)=C55(k)*(2*Ave13+W1001(k)+3*W1002(k)*Y) - 2*C44(k)*strain_inel_13(:,:,k);
        Sig12(:,:,k)=C66(k)*(2*Ave12+W1010(k)+3*W1020(k)*X) - 2*C44(k)*strain_inel_12(:,:,k);

        u1temp=cat(2, u1temp, disu1(:,:,k));
        u2temp=cat(2, u2temp, disu2(:,:,k));
        u3temp=cat(2, u3temp, disu3(:,:,k));
        Sigm11=cat(2, Sigm11, Sig11(:,:,k));
        Sigm22=cat(2, Sigm22, Sig22(:,:,k));
        Sigm33=cat(2, Sigm33, Sig33(:,:,k));
        Sigm23=cat(2, Sigm23, Sig23(:,:,k));
        Sigm13=cat(2, Sigm13, Sig13(:,:,k));
        Sigm12=cat(2, Sigm12, Sig12(:,:,k));
        Str_p_eff=cat(2, Str_p_eff, strain_inel_eff(:,:,k));
    end
    u1=cat(1, u1, u1temp);
    u2=cat(1, u2, u2temp);
    u3=cat(1, u3, u3temp);
    Sigma11=cat(1, Sigma11, Sigm11);
    Sigma22=cat(1, Sigma22, Sigm22);
    Sigma33=cat(1, Sigma33, Sigm33);
    Sigma23=cat(1, Sigma23, Sigm23);
    Sigma13=cat(1, Sigma13, Sigm13);
    Sigma12=cat(1, Sigma12, Sigm12);
    Strain_p_eff=cat(1, Strain_p_eff, Str_p_eff);
end

Sigma_HS = (Sigma11+Sigma22+Sigma33)/3;
Sigma_eff = sqrt( ((Sigma11-Sigma22).^2+(Sigma22-Sigma33).^2+(Sigma33-Sigma11).^2)/2 + 3*(Sigma12.^2+Sigma13.^2+Sigma23.^2) );
