function [eps_location, Av_e] = total_strains(m, n, h, l, G1_1, G2_1, G3_1, G1_2, G2_2, G3_2, C22, C33, C44, C55, C66, C22bar, C33bar, C55bar, C66bar, Avu, Avunext, Ave11, Ave22, Ave33, Ave23, Ave13, Ave12, ul, ur, ub, ut, array_x)

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


% %%_________________________________________________________________________
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

%%_________________________________________________________________________
Av_e11(1:NS,1)=Ave11;
Av_e22=Ave22+W2010;
Av_e33=Ave33+W3001;
Av_e23=Ave23+(W2001+W3010)/2;
Av_e13=Ave13+W1001/2;
Av_e12=Ave12+W1010/2;

Av_e=[Av_e11 Av_e22 Av_e33 2*Av_e23 2*Av_e13 2*Av_e12];
%%_________________________________________________________________________

t=length(array_x);
for i=1:t
    array_x2(:,i)=array_x(i)*h/2;
    array_x3(:,i)=array_x(i)*l/2;
end

for i=1:t
    for j=1:t
    eps_location(i,j).dir11=Av_e11;
    eps_location(i,j).dir22=Ave22+W2010+3*array_x2(:,j).*W2020;
    eps_location(i,j).dir33=Ave33+W3001+3*array_x3(:,i).*W3002;
    eps_location(i,j).dir23=Ave23+(W2001+3*array_x3(:,i).*W2002+W3010+3*array_x2(:,j).*W3020)/2;    
    eps_location(i,j).dir13=Ave13+(W1001+3*array_x3(:,i).*W1002)/2;
    eps_location(i,j).dir12=Ave12+(W1010+3*array_x2(:,j).*W1020)/2;
    end
end


% x2_1=-0.774596669241483*h/2;
% x2_2=0;
% x2_3=0.774596669241483*h/2;
% x3_1=-0.774596669241483*l/2;
% x3_2=0;
% x3_3=0.774596669241483*l/2;
%
% %%_________________________________________________________________________
% eps_11(1:NS,1)=Ave11;
%
% eps_22_1=Ave22+W2010+3*x2_1.*W2020;
% eps_22_2=Ave22+W2010+3*x2_2.*W2020;
% eps_22_3=Ave22+W2010+3*x2_3.*W2020;
%
% eps_33_1=Ave33+W3001+3*x3_1.*W3002;
% eps_33_2=Ave33+W3001+3*x3_2.*W3002;
% eps_33_3=Ave33+W3001+3*x3_3.*W3002;
%
% eps_23_11=Ave23+(W2001+3*x3_1.*W2002+W3010+3*x2_1.*W3020)/2;
% eps_23_21=Ave23+(W2001+3*x3_1.*W2002+W3010+3*x2_2.*W3020)/2;
% eps_23_31=Ave23+(W2001+3*x3_1.*W2002+W3010+3*x2_3.*W3020)/2;
% eps_23_12=Ave23+(W2001+3*x3_2.*W2002+W3010+3*x2_1.*W3020)/2;
% eps_23_22=Ave23+(W2001+3*x3_2.*W2002+W3010+3*x2_2.*W3020)/2;
% eps_23_32=Ave23+(W2001+3*x3_2.*W2002+W3010+3*x2_3.*W3020)/2;
% eps_23_13=Ave23+(W2001+3*x3_3.*W2002+W3010+3*x2_1.*W3020)/2;
% eps_23_23=Ave23+(W2001+3*x3_3.*W2002+W3010+3*x2_2.*W3020)/2;
% eps_23_33=Ave23+(W2001+3*x3_3.*W2002+W3010+3*x2_3.*W3020)/2;
%
% eps_13_1=Ave13+(W1001+3*x3_1.*W1002)/2;
% eps_13_2=Ave13+(W1001+3*x3_2.*W1002)/2;
% eps_13_3=Ave13+(W1001+3*x3_3.*W1002)/2;
%
% eps_12_1=Ave12+(W1010+3*x2_1.*W1020)/2;
% eps_12_2=Ave12+(W1010+3*x2_2.*W1020)/2;
% eps_12_3=Ave12+(W1010+3*x2_3.*W1020)/2;
% %%_________________________________________________________________________
%
%
% %%_________________________________________________________________________
% eps_location(1,1).dir11=eps_11;
% eps_location(1,2).dir11=eps_11;
% eps_location(1,3).dir11=eps_11;
% eps_location(2,1).dir11=eps_11;
% eps_location(2,2).dir11=eps_11;
% eps_location(2,3).dir11=eps_11;
% eps_location(3,1).dir11=eps_11;
% eps_location(3,2).dir11=eps_11;
% eps_location(3,3).dir11=eps_11;
%
% eps_location(1,1).dir22=eps_22_1;
% eps_location(1,2).dir22=eps_22_2;
% eps_location(1,3).dir22=eps_22_3;
% eps_location(2,1).dir22=eps_22_1;
% eps_location(2,2).dir22=eps_22_2;
% eps_location(2,3).dir22=eps_22_3;
% eps_location(3,1).dir22=eps_22_1;
% eps_location(3,2).dir22=eps_22_2;
% eps_location(3,3).dir22=eps_22_3;
%
% eps_location(1,1).dir33=eps_33_1;
% eps_location(1,2).dir33=eps_33_1;
% eps_location(1,3).dir33=eps_33_1;
% eps_location(2,1).dir33=eps_33_2;
% eps_location(2,2).dir33=eps_33_2;
% eps_location(2,3).dir33=eps_33_2;
% eps_location(3,1).dir33=eps_33_3;
% eps_location(3,2).dir33=eps_33_3;
% eps_location(3,3).dir33=eps_33_3;
%
% eps_location(1,1).dir23=eps_23_11;
% eps_location(1,2).dir23=eps_23_21;
% eps_location(1,3).dir23=eps_23_31;
% eps_location(2,1).dir23=eps_23_12;
% eps_location(2,2).dir23=eps_23_22;
% eps_location(2,3).dir23=eps_23_32;
% eps_location(3,1).dir23=eps_23_13;
% eps_location(3,2).dir23=eps_23_23;
% eps_location(3,3).dir23=eps_23_33;
%
% eps_location(1,1).dir13=eps_13_1;
% eps_location(1,2).dir13=eps_13_1;
% eps_location(1,3).dir13=eps_13_1;
% eps_location(2,1).dir13=eps_13_2;
% eps_location(2,2).dir13=eps_13_2;
% eps_location(2,3).dir13=eps_13_2;
% eps_location(3,1).dir13=eps_13_3;
% eps_location(3,2).dir13=eps_13_3;
% eps_location(3,3).dir13=eps_13_3;
%
% eps_location(1,1).dir12=eps_12_1;
% eps_location(1,2).dir12=eps_12_2;
% eps_location(1,3).dir12=eps_12_3;
% eps_location(2,1).dir12=eps_12_1;
% eps_location(2,2).dir12=eps_12_2;
% eps_location(2,3).dir12=eps_12_3;
% eps_location(3,1).dir12=eps_12_1;
% eps_location(3,2).dir12=eps_12_2;
% eps_location(3,3).dir12=eps_12_3;
% %%_________________________________________________________________________
