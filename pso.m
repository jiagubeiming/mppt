function D = PSO(Vpv,Ipv)

%变量
persistent u;%粒子
persistent dcurrent;%当前占空比
persistent pbest;%个体
persistent p;%功率
persistent dc;%占空比
persistent v;%速度
persistent counter;%计数
persistent gbest;%群体
%初始化
if(isempty(counter))
    counter=0;
end
if(isempty(dcurrent))
    dcurrent=0.5;
end
if(isempty(gbest))
    gbest=0.5;
end
if(isempty(p))
    p=zeros(4,1);
end

if(isempty(v))
    v=zeros(4,1);
end

if(isempty(pbest))
    pbest=zeros(4,1);
end
if(isempty(u))
    u=0;
end
if(isempty(dc))%[0;0.3;0.5;0.9]
    dc=zeros(4,1);
    dc(1)=0;
    dc(2)=0.3;
    dc(3)=0.5;
    dc(4)=0.9;
end

%延时
if(counter>=1 && counter<300)%
    D=dcurrent;
    counter=counter+1;
    return;%跳出循环
end
counter=0;
%计算每一个粒子的功率值
%将结果与前一个比较，若优则更新
if(u>=1 && u<=4)
    if((Vpv*Ipv)>p(u))
        p(u)=Vpv*Ipv;
        pbest(u)=dcurrent;
    end
end
u=u+1;

if(u==6)
    u=1;
end
if(u>=1 && u<=4)
    D=dc(u);
    dcurrent=D;
    counter=1;
    return;
elseif(u==5 )
    [m,i]=max(p);%找出最大值
    gbest=pbest(i);%最大功率粒子的位置（占空比）
    D=gbest;
    dcurrent=D;
    counter=1;

    %%更新速度
    v(1)=updatevelocity(v(1),pbest(1),dc(1),gbest)
    v(2)=updatevelocity(v(2),pbest(2),dc(2),gbest)
    v(3)=updatevelocity(v(3),pbest(3),dc(3),gbest)
    v(4)=updatevelocity(v(4),pbest(4),dc(4),gbest)
    %更新占空比
    dc(1)=updateduty(dc(1),v(1))
    dc(2)=updateduty(dc(2),v(2))
    dc(3)=updateduty(dc(3),v(3))
    dc(4)=updateduty(dc(4),v(4))
    return;
else
    D=0.1 
end
end
%
function vfinal=updatevelocity(velocity,pobest,d,gwbest)
% 
w=0.4;          % 惯性权重 %w=(1.1-(pobest/gwbest));
c1=1.2;          % 自我学习因子
c2=2;          % 群体学习因子
vfinal = (w*velocity)+(c1*rand(1)*(pobest-d))+(c2*rand(1)*(gwbest-d));
end

function dfinal=updateduty(d,velocity)
dup=d+velocity;
if(dup>1)
    dfinal=1;
elseif(dup<0)
    dfinal=0;
else
    dfinal=dup;
end
end


% 使用收缩系数
% % Constriction Coefficients
% phi1=2.05;
% phi2=2.05;
% phi=phi1+phi2;
% chi=2/(phi-2+sqrt(phi^2-4*phi));
% w=chi;          % 惯性权重
% wdamp=1;        % 惯性权重阻尼比
% c1=chi*phi1;    % 自我学习因子
% c2=chi*phi2;    % 群体学习因子