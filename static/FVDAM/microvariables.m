function [W2010, W3010, W2001, W3001, W1010, W1001] = microvariables(m, n, h, l, Avu, Avunext, ul, ur, ub, ut);

% Calculate the microvariables for each subcell using the surface averaged displacements

%%_________________________________________________________________________
NS=m*n;
u221=zeros(1,NS);
u222=zeros(1,NS);
u231=zeros(1,NS);
u232=zeros(1,NS);
u321=zeros(1,NS);
u322=zeros(1,NS);
u331=zeros(1,NS);
u332=zeros(1,NS);
u211=zeros(1,NS);
u212=zeros(1,NS);
u311=zeros(1,NS);
u312=zeros(1,NS);

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
W1010=(u212-u211)./h;
W1001=(u312-u311)./l;
W2010=(u222-u221)./h;
W3010=(u232-u231)./h;
W2001=(u322-u321)./l;
W3001=(u332-u331)./l;
% for i=1:m*n
%     W1010(i)=(u212(i)-u211(i))/h(i);
%     W1001(i)=(u312(i)-u311(i))/l(i);
%     W2010(i)=(u222(i)-u221(i))/h(i);
%     W3010(i)=(u232(i)-u231(i))/h(i);
%     W2001(i)=(u322(i)-u321(i))/l(i);
%     W3001(i)=(u332(i)-u331(i))/l(i);
% end
%%_________________________________________________________________________
