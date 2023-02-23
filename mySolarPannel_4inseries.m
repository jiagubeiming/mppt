clear
clc
 
%-----------------------------------------------
%-----------------------------------------------
%pannel in series
%first pannel
S_1=1000;
Tair_1=25;
 
Sref=1000;  %1000W/m^2
Tref=25;    %25degree celcius
 
Uoc=36.8;
Um=30;
Isc=8.83;
Im=8.3;
 
a=0.00255;
b=0.55;
c=0.00285;
 
T_1 = Tair_1 + 0.028*S_1;
T_delta_1 = T_1 - Tref;
S_delta_1 = S_1/Sref - 1;
 
Isc_comp_1 = Isc*S_1/Sref*(1+a*T_delta_1);
Uoc_comp_1 = Uoc*(1-c*T_delta_1)*log(exp(1)+b*S_delta_1);
Im_comp_1  = Im*S_1/Sref*(1+a*T_delta_1);
Um_comp_1  = Um*(1-c*T_delta_1)*log(exp(1)+b*S_delta_1);
 
C2_1=(Um_comp_1/Uoc_comp_1-1)*(log(1-Im_comp_1/Isc_comp_1))^(-1);
C1_1=(1-Im_comp_1/Isc_comp_1)*exp(-Um_comp_1/(C2_1*Uoc_comp_1));
 
U_1=0:0.01:Uoc_comp_1;
Iph_1=Isc_comp_1*(1-C1_1*(exp(U_1/(C2_1*Uoc_comp_1))-1));
 
 
%-----------------------------------------------
%second pannel
S_2=1000;
Tair_2=25;
 
Sref=1000;  %1000W/m^2
Tref=25;    %25degree celcius
 
Uoc=36.8;
Um=30;
Isc=8.83;
Im=8.3;
 
a=0.00255;
b=0.55;
c=0.00285;
 
T_2 = Tair_2 + 0.028*S_2;
T_delta_2 = T_2 - Tref;
S_delta_2 = S_2/Sref - 1;
 
Isc_comp_2 = Isc*S_2/Sref*(1+a*T_delta_2);
Uoc_comp_2 = Uoc*(1-c*T_delta_2)*log(exp(1)+b*S_delta_2);
Im_comp_2  = Im*S_2/Sref*(1+a*T_delta_2);
Um_comp_2  = Um*(1-c*T_delta_2)*log(exp(1)+b*S_delta_2);
 
C2_2=(Um_comp_2/Uoc_comp_2-1)*(log(1-Im_comp_2/Isc_comp_2))^(-1);
C1_2=(1-Im_comp_2/Isc_comp_2)*exp(-Um_comp_2/(C2_2*Uoc_comp_2));
 
U_2=0:0.01:Uoc_comp_2;
Iph_2=Isc_comp_2*(1-C1_2*(exp(U_2/(C2_2*Uoc_comp_2))-1));
 
 
 
%-----------------------------------------------
%third pannel
S_3=800;
Tair_3=25;
 
Sref=1000;  %1000W/m^2
Tref=25;    %25degree celcius
 
Uoc=36.8;
Um=30;
Isc=8.83;
Im=8.3;
 
a=0.00255;
b=0.55;
c=0.00285;
 
T_3 = Tair_3 + 0.028*S_3;
T_delta_3 = T_3 - Tref;
S_delta_3 = S_3/Sref - 1;
 
Isc_comp_3 = Isc*S_3/Sref*(1+a*T_delta_3);
Uoc_comp_3 = Uoc*(1-c*T_delta_3)*log(exp(1)+b*S_delta_3);
Im_comp_3  = Im*S_3/Sref*(1+a*T_delta_3);
Um_comp_3  = Um*(1-c*T_delta_3)*log(exp(1)+b*S_delta_3);
 
C2_3=(Um_comp_3/Uoc_comp_3-1)*(log(1-Im_comp_3/Isc_comp_3))^(-1);
C1_3=(1-Im_comp_3/Isc_comp_3)*exp(-Um_comp_3/(C2_3*Uoc_comp_3));
 
U_3=0:0.01:Uoc_comp_3;
Iph_3=Isc_comp_3*(1-C1_3*(exp(U_3/(C2_3*Uoc_comp_3))-1));
 
 
 
%-----------------------------------------------
%forth pannel
S_4=500;
Tair_4=25;
 
Sref=1000;  %1000W/m^2
Tref=25;    %25degree celcius
 
Uoc=36.8;
Um=30;
Isc=8.83;
Im=8.3;
 
a=0.00255;
b=0.55;
c=0.00285;
 
T_4 = Tair_4 + 0.028*S_4;
T_delta_4 = T_4 - Tref;
S_delta_4 = S_4/Sref - 1;
 
Isc_comp_4 = Isc*S_4/Sref*(1+a*T_delta_4);
Uoc_comp_4 = Uoc*(1-c*T_delta_4)*log(exp(1)+b*S_delta_4);
Im_comp_4  = Im*S_4/Sref*(1+a*T_delta_4);
Um_comp_4  = Um*(1-c*T_delta_4)*log(exp(1)+b*S_delta_4);
 
C2_4=(Um_comp_4/Uoc_comp_4-1)*(log(1-Im_comp_4/Isc_comp_4))^(-1);
C1_4=(1-Im_comp_4/Isc_comp_4)*exp(-Um_comp_4/(C2_4*Uoc_comp_4));
 
U_4=0:0.01:Uoc_comp_4;
Iph_4=Isc_comp_4*(1-C1_4*(exp(U_4/(C2_4*Uoc_comp_4))-1));
 
%{
plot(U_1,Iph_1)
hold on
plot(U_2,Iph_2)
hold on
plot(U_3,Iph_3)
hold on
plot(U_4,Iph_4)
%}
%-----------------------------------------------
% 4 in series
% U=C2*Uoc*log((Isc-I)/(C1*Isc)+1)
%Iph_1 > Iph_2 > Iph_3
 
 
U_s = 0:0.01:Uoc_comp_1*4;
Iph_s=Isc_comp_1*(1-C1_1*(exp(U_s/(C2_1*Uoc_comp_1*4))-1));
plot(U_s,Iph_s,'k')
 
U_ss = zeros(size(U_1)+size(U_2)+size(U_3)+size(U_4));
Iph_ss = U_ss;
%for i = 1 : size(U_1)(2)
%    U_ss(i) = U_1(i);
%    Iph_ss(i) = Iph_1(i);
%end
 
 
for i = 1 : length(U_1)
    if Iph_1(i)>=Iph_2(1)
        U_ss(i) = U_1(i);
        Iph_ss(i) = Iph_1(i);
        step1 = i;
    else
        break;
    end
end
for i = 1 : length(U_2)
    if Iph_2(i)>Iph_3(1)
        U_ss(step1+i) = U_2(i) + U_ss(step1);
        Iph_ss(step1+i) = Iph_2(i);
        step2 = step1+i;
    else
        break;
    end
end
for i = 1 : length(U_3)
    if Iph_3(i)>Iph_4(1)
        U_ss(step2+i) = U_3(i) + U_ss(step2);
        Iph_ss(step2+i) = Iph_3(i);
        step3 = step2+i;
    else
        break;
    end
end
for i = 1 : length(U_4)
    U_ss(step3+i) = U_ss(step3) + U_4(i);
    Iph_ss(step3+i) = Iph_4(i);
end
 
%plot(U_1,Iph_1)               I-U
%plot(U_2,Iph_2)
%plot(U_ss,Iph_ss,'+')
%plot(U_ss,Iph_ss,'+')
 
figure(1)
plot(U_ss,Iph_ss,'+')
xlabel('U/V')
ylabel('I/A')
title('U-I')
figure(2)
P_ss = U_ss .* Iph_ss;
plot(U_ss,P_ss)
xlabel('U/V')
ylabel('P/W')
title('U-W')
hold on
P_1 = U_1*4 .* Iph_1;
plot(U_1*4,P_1)