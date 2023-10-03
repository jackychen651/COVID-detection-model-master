clc
clear
close all;
h1=0.001;v=2000000;N=25000000;
h2=0.15;
z1=0.001;z2=0.001;
p=zeros(60,1);I=zeros(60,1);S=zeros(60,1);a=zeros(60,1);t=zeros(60,1);
k=zeros(60,1);J=zeros(60,1);
k_spr1=zeros(60,1);k_spr2=zeros(60,1);%传染率
k_spr1(1)=0.08;k_spr2(1)=0.001;%
a(1)=0.3;
I(1) = 0.0003;
S(1) = 0.9996;p(1)=(1-S(1))*20;
J(1)=1-I(1)-S(1);
i=0;i_f=0;
while (S(i+1)<0.999996) && (i<59)
    i=i+1;
    %局部变量复制（同一个周期）
    p_local=p(i);S_local=S(i);I_local=I(i);a_local=a(i);
    k_sprloc1=k_spr1(i);k_sprloc2=k_spr2(i);
    %得出该阶段的最佳混检率
    fun = @(x)1.5+1/x-1.5*(1-p_local)^x;
    x0 = 3;
    [k_local1,eval] = fminunc(fun,x0);
    k_local=round(k_local1);k(i)=k_local;
    M=N*(1.5+1/k_local-1.5*(1-p_local)^(k_local));t(i)=ceil(M/v);
    a(i+1)=2*(1-p_local)^k_local*N/(v*t(i))*a(1);
    %ode求解
    t_loc=t(i);
    tspan = [0 t_loc];
    y0 = [S_local I_local];
    [t1, y] = ode45(@(t,y)odefun(t,y,k_sprloc1,k_sprloc2,h1,h2,a_local), tspan, y0);
    %更新参数
    Step=max(0,y(end,1));S(i+1)=min(1,Step);
    Itep=max(0,y(end,2));I(i+1)=min(1-Step,Itep);
    J(i+1) = max(0,1-S(i+1)-I(i+1));
    p(i+1)=min(0.012-(i/6000),(1-S(i+1))*20);
    %k_spr1(i+1)=max(0,(y(end,2)-y(end-1,2))*6/y(end-1,2)+h1);
    %k_spr2(i+1)=max(0,(y(end-1,1)+y(end-1,2)-y(end,1)-y(end,2))*2/(1-y(end-1,1)-y(end-1,2)+h2));
    k_spr1(i+1)=k_spr1(i);k_spr2(i+1)=k_spr2(i);
    %clear t1;clear y;
end
 T=zeros(60,1);
for j=1:60
    T(j)=j;
end
I_f=zeros(60,3);S_f=zeros(60,3);t_f=zeros(60,3);k_f=zeros(60,3);J_f=zeros(60,3);
for tot=1:i
    I_f(tot,1)=I(tot); S_f(tot,1)=S(tot); J_f(tot,1)=J(tot);
    t_f(tot,1)=t(tot); k_f(tot,1)=k(tot);
end
i_f=i;
%第二轮
i=0;p(1)=(1-S(1))*20*2;
while (S(i+1)<0.999996) && (i<59)
    i=i+1;
    %局部变量复制（同一个周期）
    p_local=p(i);S_local=S(i);I_local=I(i);a_local=a(i);
    k_sprloc1=k_spr1(i);k_sprloc2=k_spr2(i);
    %得出该阶段的最佳混检率
    fun = @(x)1.5+1/x-1.5*(1-p_local)^x;
    x0 = 3;
    [k_local1,eval] = fminunc(fun,x0);
    k_local=round(k_local1);k(i)=k_local;
    M=N*(1.5+1/k_local-1.5*(1-p_local)^(k_local));t(i)=ceil(M/v);
    a(i+1)=2*(1-p_local)^k_local*N/(v*t(i))*a(1);
    %ode求解
    t_loc=t(i);
    tspan = [0 t_loc];
    y0 = [S_local I_local];
    [t1, y] = ode45(@(t,y)odefun(t,y,k_sprloc1,k_sprloc2,h1,h2,a_local), tspan, y0);
    %更新参数
    Step=max(0,y(end,1));S(i+1)=min(1,Step);
    Itep=max(0,y(end,2));I(i+1)=min(1-Step,Itep);
    J(i+1) = max(0,1-S(i+1)-I(i+1));
    p(i+1)=min(0.012-(i/6000),(1-S(i+1))*20);
    %k_spr1(i+1)=max(0,(y(end,2)-y(end-1,2))*6/y(end-1,2)+h1);
    %k_spr2(i+1)=max(0,(y(end-1,1)+y(end-1,2)-y(end,1)-y(end,2))*2/(1-y(end-1,1)-y(end-1,2)+h2));
    k_spr1(i+1)=k_spr1(i);k_spr2(i+1)=k_spr2(i);
    %clear t1;clear y;
end
for j=1:60
    T(j)=j;
end
for tot=1:i
    I_f(tot,2)=I(tot); S_f(tot,2)=S(tot); J_f(tot,2)=J(tot);
    t_f(tot,2)=t(tot); k_f(tot,2)=k(tot);
end
if i_f>i
    i_f=i;
end
%第三轮
i=0;p(1)=(1-S(1))*20*0.5;
while (S(i+1)<0.999996) && (i<59)
    i=i+1;
    %局部变量复制（同一个周期）
    p_local=p(i);S_local=S(i);I_local=I(i);a_local=a(i);
    k_sprloc1=k_spr1(i);k_sprloc2=k_spr2(i);
    %得出该阶段的最佳混检率
    fun = @(x)1.5+1/x-1.5*(1-p_local)^x;
    x0 = 3;
    [k_local1,eval] = fminunc(fun,x0);
    k_local=round(k_local1);k(i)=k_local;
    M=N*(1.5+1/k_local-1.5*(1-p_local)^(k_local));t(i)=ceil(M/v);
    a(i+1)=2*(1-p_local)^k_local*N/(v*t(i))*a(1);
    %ode求解
    t_loc=t(i);
    tspan = [0 t_loc];
    y0 = [S_local I_local];
    [t1, y] = ode45(@(t,y)odefun(t,y,k_sprloc1,k_sprloc2,h1,h2,a_local), tspan, y0);
    %更新参数
    Step=max(0,y(end,1));S(i+1)=min(1,Step);
    Itep=max(0,y(end,2));I(i+1)=min(1-Step,Itep);
    J(i+1) = max(0,1-S(i+1)-I(i+1));
    p(i+1)=min(0.012-(i/6000),(1-S(i+1))*20);
    %k_spr1(i+1)=max(0,(y(end,2)-y(end-1,2))*6/y(end-1,2)+h1);
    %k_spr2(i+1)=max(0,(y(end-1,1)+y(end-1,2)-y(end,1)-y(end,2))*2/(1-y(end-1,1)-y(end-1,2)+h2));
    k_spr1(i+1)=k_spr1(i);k_spr2(i+1)=k_spr2(i);
    %clear t1;clear y;
end
for j=1:60
    T(j)=j;
end
for tot=1:i
    I_f(tot,3)=I(tot); S_f(tot,3)=S(tot); J_f(tot,3)=J(tot);
    t_f(tot,3)=t(tot); k_f(tot,3)=k(tot);
end
if i_f>i
    i_f=i;
end

%画图
subplot(2,2,1);
plot(T(1:i_f),S_f(1:i_f,1),'-o',T(1:i_f),S_f(1:i_f,2),'-o',T(1:i_f),S_f(1:i_f,3),'-o');
legend('Norm:S(t)','Less:S1(t)','More:S2(t)'); 
ylabel('Proportion of Population ');
xlabel('Period T');
hold on;
subplot(2,2,2)
plot(T(1:i),I_f(1:i,1),T(1:i),I_f(1:i,2),'g',T(1:i),I_f(1:i,3));
legend('Norm :I(t)','Less :I1(t)','More:I2(t)','Best'); 
ylabel('Proportion of Population ');
xlabel('Period T');
hold on;
subplot(2,2,3)
plot(T(1:i),J_f(1:i,1),T(1:i),J_f(1:i,2),'g',T(1:i),J_f(1:i,3));
legend('Norm :J(t)','Less :J1(t)','More:J2(t)','Best'); 
ylabel('Proportion of Population ');
xlabel('Period T');
hold on;

hold on;
subplot(2,2,4)
plot(T(1:i),k_f(1:i,1),'-.',T(1:i),k_f(1:i,2),'-.',T(1:i),k_f(1:i,3),'-.');
legend('Norm :K(t)','Less :K1(t)','More:K2(t)','Best'); 
ylabel('Number');
xlabel('Period T');







%y(1)是S,y（2）是I,y(3)是J，实际仅需两个方程即可;z1是死亡率
function dydt = odefun(t,y,k1,k2,h1,h2,a)
dydt = zeros(2,1);
dydt(1) = -(k1*y(2)+(1-y(1)-y(2))*k2)*y(1)+h1*y(2)+h2*(1-y(1)-y(2));
dydt(2) = (k1*y(2)+k2*(1-y(1)-y(2)))*y(1)-(a+h1)*y(2);
end 
% function dydt = odefun(t,y,k1,k2,a)
% dydt = zeros(2,1);
% dydt(1) = -(k1*y(2)+(1-y(1)-y(2)))*k2*y(1);
% dydt(2) =(k1*y(2)+k2*(1-y(1)-y(2)))*y(1)-a*y(2);
% end 