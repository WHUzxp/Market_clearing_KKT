%�����ǵ綯�����ɵ���Ǳ���ĳ��վ���׶��г�Ͷ����ԡ�һ����س�����Ϊ�������������������
%�������ܣ�ʵʱ�ɵ���Ǳ�����㣬��ǰ�ɵ���Ǳ�����㣬���վ��ǰͶ����ԣ����վʵʱͶ����ԣ����������ȷ���
%����ѧ������ʹ�ã���ע���������Ȩ
%authors:X.P. Zhan. GitHub:WHU_zxp
%���л�����MATLAB R2014a,��ҪYALMIP�������GUROBI�����
clear
clc
close all
%% ��һ��(�������ݻ�ȡ���������ؿ���������ɵ���Ǳ�������Լ����������)
%%%���ؿ������%%%
%I�೵(800)��ҹ���磬˽�ҳ���II�೵(400)��ҹ���磬��Լ����III�೵(800)�������磬�ϰ��塣
rng(1);%�̶����������(α�������֤���ؿ�������Ľ���ȶ�)
Ta_EV1=normrnd(18,2,[1001,1000]);%�綯����1��ͣ��ʱ�����
rng(1);Tl_EV1=normrnd(8,2,[1001,1000]);%�綯����1���뿪ʱ�����
rng(1);S0_EV1=unifrnd(0.4,0.6,[1001,1000]);%�綯����1�ĳ�ʼSOC����
rng(1);Ta_EV2=normrnd(21,1,[1001,500]);
rng(1);Tl_EV2=normrnd(7,1,[1001,500]);
rng(1);S0_EV2=unifrnd(0.2,0.4,[1001,500]);
rng(1);Ta_EV3=normrnd(9,sqrt(2),[1001,1000]);
rng(1);Tl_EV3=normrnd(17,sqrt(2),[1001,1000]);
rng(1);S0_EV3=unifrnd(0.4,0.6,[1001,1000]);
rng(1);EV1_CS1=round(unifrnd(180,220,[1001,1]));%���վ1��I�೵������
rng(1);EV1_CS2=round(unifrnd(180,220,[1001,1]));%���վ2��I�೵������
rng(1);EV1_CS3=zeros(1001,1);%���վ3��I�೵������
rng(1);EV1_CS4=round(unifrnd(380,420,[1001,1]));%���վ4��I�೵������
rng(1);EV2_CS1=round(unifrnd(190,210,[1001,1]));
rng(1);EV2_CS2=round(unifrnd(80,120,[1001,1]));
rng(1);EV2_CS3=round(unifrnd(90,110,[1001,1]));
rng(1);EV2_CS4=zeros(1001,1);
rng(1);EV3_CS1=zeros(1001,1);
rng(1);EV3_CS2=round(unifrnd(380,420,[1001,1]));
rng(1);EV3_CS3=round(unifrnd(380,420,[1001,1]));
rng(1);EV3_CS4=zeros(1001,1);
%������������08:00��ʼ���㣬������Χ�Ľ��н�β(��β��̬�ֲ�)
Ta_EV1=Ta_EV1-8;Tl_EV1=Tl_EV1+24-8;Ta_EV2=Ta_EV2-8;Tl_EV2=Tl_EV2+24-8;Ta_EV3=Ta_EV3-8;Tl_EV3=Tl_EV3-8;S0_EV1=32*S0_EV1;S0_EV2=32*S0_EV2;S0_EV3=32*S0_EV3;
DA_Ta_EV1=ceil(Ta_EV1(1:1000,:));DA_Ta_EV2=ceil(Ta_EV2(1:1000,:));DA_Ta_EV3=ceil(Ta_EV3(1:1000,:));%��ǰ��ɢ��
DA_Tl_EV1=floor(Tl_EV1(1:1000,:));DA_Tl_EV2=floor(Tl_EV2(1:1000,:));DA_Tl_EV3=floor(Tl_EV3(1:1000,:));
DA_S0_EV1=S0_EV1(1:1000,:);DA_S0_EV2=S0_EV2(1:1000,:);DA_S0_EV3=S0_EV3(1:1000,:);
RT_Ta_EV1=ceil(4*Ta_EV1(1001,:));RT_Ta_EV2=ceil(4*Ta_EV2(1001,:));RT_Ta_EV3=ceil(4*Ta_EV3(1001,:));%ʵʱ��ɢ��
RT_Tl_EV1=floor(4*Tl_EV1(1001,:));RT_Tl_EV2=floor(4*Tl_EV2(1001,:));RT_Tl_EV3=floor(4*Tl_EV3(1001,:));
RT_S0_EV1=S0_EV1(1001,:);RT_S0_EV2=S0_EV2(1001,:);RT_S0_EV3=S0_EV3(1001,:);
DA_Ta_EV1(find(DA_Ta_EV1<1))=1;DA_Ta_EV2(find(DA_Ta_EV2<1))=1;DA_Ta_EV3(find(DA_Ta_EV3<1))=1;%Ta����
RT_Ta_EV1(find(RT_Ta_EV1<1))=1;RT_Ta_EV2(find(RT_Ta_EV2<1))=1;RT_Ta_EV3(find(RT_Ta_EV3<1))=1;
DA_Ta_EV1(find(DA_Ta_EV1>24))=24;DA_Ta_EV2(find(DA_Ta_EV2>24))=24;DA_Ta_EV3(find(DA_Ta_EV3>24))=24;%Ta��β
RT_Ta_EV1(find(RT_Ta_EV1>96))=96;RT_Ta_EV2(find(RT_Ta_EV2>96))=96;RT_Ta_EV3(find(RT_Ta_EV3>96))=96;
DA_Tl_EV1(find(DA_Tl_EV1<1))=1;DA_Tl_EV2(find(DA_Tl_EV2<1))=1;DA_Tl_EV3(find(DA_Tl_EV3<1))=1;%Tl����
RT_Tl_EV1(find(RT_Tl_EV1<1))=1;RT_Tl_EV2(find(RT_Tl_EV2<1))=1;RT_Tl_EV3(find(RT_Tl_EV3<1))=1;
DA_Tl_EV1(find(DA_Tl_EV1>24))=24;DA_Tl_EV2(find(DA_Tl_EV2>24))=24;DA_Tl_EV3(find(DA_Tl_EV3>24))=24;%Tl��β
RT_Tl_EV1(find(RT_Tl_EV1>96))=96;RT_Tl_EV2(find(RT_Tl_EV2>96))=96;RT_Tl_EV3(find(RT_Tl_EV3>96))=96;
index=find(6.6*(DA_Tl_EV1-DA_Ta_EV1+1)<32*0.9-DA_S0_EV1);%��ǰ������������������������ֵ
for i=1:length(index)
    DA_Ta_EV1(index(i))=10;DA_Tl_EV1(index(i))=24;
end
index=find(6.6*(DA_Tl_EV2-DA_Ta_EV2+1)<32*0.9-DA_S0_EV2);
for i=1:length(index)
    DA_Ta_EV2(index(i))=13;DA_Tl_EV2(index(i))=23;
end
index=find(6.6*(DA_Tl_EV3-DA_Ta_EV3+1)<32*0.9-DA_S0_EV3);
for i=1:length(index)
    DA_Ta_EV3(index(i))=1;DA_Tl_EV3(index(i))=9;
end
index=find(0.25*6.6*(RT_Tl_EV1-RT_Ta_EV1+1)<32*0.9-RT_S0_EV1);%ʵʱ������������������������ֵ
for i=1:length(index)
    RT_Ta_EV1(index(i))=40;RT_Tl_EV1(index(i))=96;
end
index=find(0.25*6.6*(RT_Tl_EV2-RT_Ta_EV2+1)<32*0.9-RT_S0_EV2);
for i=1:length(index)
    RT_Ta_EV2(index(i))=52;RT_Tl_EV2(index(i))=92;
end
index=find(0.25*6.6*(RT_Tl_EV3-RT_Ta_EV3+1)<32*0.9-RT_S0_EV3);
for i=1:length(index)
    RT_Ta_EV3(index(i))=4;RT_Tl_EV3(index(i))=36;
end
%%%��������%%%
for i=1:1000
    data_CS1(i).Ta=[DA_Ta_EV1(i,1:EV1_CS1(i)),DA_Ta_EV2(i,1:EV2_CS1(i)),DA_Ta_EV3(i,1:EV3_CS1(i))];
    data_CS2(i).Ta=[DA_Ta_EV1(i,1+EV1_CS1(i):EV1_CS1(i)+EV1_CS2(i)),DA_Ta_EV2(i,1+EV2_CS1(i):EV2_CS1(i)+EV2_CS2(i)),DA_Ta_EV3(i,1+EV3_CS1(i):EV3_CS1(i)+EV3_CS2(i))];
    data_CS3(i).Ta=[DA_Ta_EV1(i,1+EV1_CS1(i)+EV1_CS2(i):EV1_CS1(i)+EV1_CS2(i)+EV1_CS3(i)),DA_Ta_EV2(i,1+EV2_CS1(i)+EV2_CS2(i):EV2_CS1(i)+EV2_CS2(i)+EV2_CS3(i)),DA_Ta_EV3(i,1+EV3_CS1(i)+EV3_CS2(i):EV3_CS1(i)+EV3_CS2(i)+EV3_CS3(i))];
    data_CS4(i).Ta=[DA_Ta_EV1(i,1+EV1_CS1(i)+EV1_CS2(i)+EV1_CS3(i):EV1_CS1(i)+EV1_CS2(i)+EV1_CS3(i)+EV1_CS4(i)),DA_Ta_EV2(i,1+EV2_CS1(i)+EV2_CS2(i)+EV2_CS3(i):EV2_CS1(i)+EV2_CS2(i)+EV2_CS3(i)+EV2_CS4(i)),DA_Ta_EV3(i,1+EV3_CS1(i)+EV3_CS2(i)+EV3_CS3(i):EV3_CS1(i)+EV3_CS2(i)+EV3_CS3(i)+EV3_CS4(i))];
    data_CS1(i).Tl=[DA_Tl_EV1(i,1:EV1_CS1(i)),DA_Tl_EV2(i,1:EV2_CS1(i)),DA_Tl_EV3(i,1:EV3_CS1(i))];
    data_CS2(i).Tl=[DA_Tl_EV1(i,1+EV1_CS1(i):EV1_CS1(i)+EV1_CS2(i)),DA_Tl_EV2(i,1+EV2_CS1(i):EV2_CS1(i)+EV2_CS2(i)),DA_Tl_EV3(i,1+EV3_CS1(i):EV3_CS1(i)+EV3_CS2(i))];
    data_CS3(i).Tl=[DA_Tl_EV1(i,1+EV1_CS1(i)+EV1_CS2(i):EV1_CS1(i)+EV1_CS2(i)+EV1_CS3(i)),DA_Tl_EV2(i,1+EV2_CS1(i)+EV2_CS2(i):EV2_CS1(i)+EV2_CS2(i)+EV2_CS3(i)),DA_Tl_EV3(i,1+EV3_CS1(i)+EV3_CS2(i):EV3_CS1(i)+EV3_CS2(i)+EV3_CS3(i))];
    data_CS4(i).Tl=[DA_Tl_EV1(i,1+EV1_CS1(i)+EV1_CS2(i)+EV1_CS3(i):EV1_CS1(i)+EV1_CS2(i)+EV1_CS3(i)+EV1_CS4(i)),DA_Tl_EV2(i,1+EV2_CS1(i)+EV2_CS2(i)+EV2_CS3(i):EV2_CS1(i)+EV2_CS2(i)+EV2_CS3(i)+EV2_CS4(i)),DA_Tl_EV3(i,1+EV3_CS1(i)+EV3_CS2(i)+EV3_CS3(i):EV3_CS1(i)+EV3_CS2(i)+EV3_CS3(i)+EV3_CS4(i))];
    data_CS1(i).S0=[DA_S0_EV1(i,1:EV1_CS1(i)),DA_S0_EV2(i,1:EV2_CS1(i)),DA_S0_EV3(i,1:EV3_CS1(i))];
    data_CS2(i).S0=[DA_S0_EV1(i,1+EV1_CS1(i):EV1_CS1(i)+EV1_CS2(i)),DA_S0_EV2(i,1+EV2_CS1(i):EV2_CS1(i)+EV2_CS2(i)),DA_S0_EV3(i,1+EV3_CS1(i):EV3_CS1(i)+EV3_CS2(i))];
    data_CS3(i).S0=[DA_S0_EV1(i,1+EV1_CS1(i)+EV1_CS2(i):EV1_CS1(i)+EV1_CS2(i)+EV1_CS3(i)),DA_S0_EV2(i,1+EV2_CS1(i)+EV2_CS2(i):EV2_CS1(i)+EV2_CS2(i)+EV2_CS3(i)),DA_S0_EV3(i,1+EV3_CS1(i)+EV3_CS2(i):EV3_CS1(i)+EV3_CS2(i)+EV3_CS3(i))];
    data_CS4(i).S0=[DA_S0_EV1(i,1+EV1_CS1(i)+EV1_CS2(i)+EV1_CS3(i):EV1_CS1(i)+EV1_CS2(i)+EV1_CS3(i)+EV1_CS4(i)),DA_S0_EV2(i,1+EV2_CS1(i)+EV2_CS2(i)+EV2_CS3(i):EV2_CS1(i)+EV2_CS2(i)+EV2_CS3(i)+EV2_CS4(i)),DA_S0_EV3(i,1+EV3_CS1(i)+EV3_CS2(i)+EV3_CS3(i):EV3_CS1(i)+EV3_CS2(i)+EV3_CS3(i)+EV3_CS4(i))];
end
data_CS1(1001).Ta=[RT_Ta_EV1(1:EV1_CS1(1001)),RT_Ta_EV2(1:EV2_CS1(1001)),RT_Ta_EV3(1:EV3_CS1(1001))];
data_CS2(1001).Ta=[RT_Ta_EV1(1+EV1_CS1(1001):EV1_CS1(1001)+EV1_CS2(1001)),RT_Ta_EV2(1+EV2_CS1(1001):EV2_CS1(1001)+EV2_CS2(1001)),RT_Ta_EV3(1+EV3_CS1(1001):EV3_CS1(1001)+EV3_CS2(1001))];
data_CS3(1001).Ta=[RT_Ta_EV1(1+EV1_CS1(1001)+EV1_CS2(1001):EV1_CS1(1001)+EV1_CS2(1001)+EV1_CS3(1001)),RT_Ta_EV2(1+EV2_CS1(1001)+EV2_CS2(1001):EV2_CS1(1001)+EV2_CS2(1001)+EV2_CS3(1001)),RT_Ta_EV3(1+EV3_CS1(1001)+EV3_CS2(1001):EV3_CS1(1001)+EV3_CS2(1001)+EV3_CS3(1001))];
data_CS4(1001).Ta=[RT_Ta_EV1(1+EV1_CS1(1001)+EV1_CS2(1001)+EV1_CS3(1001):EV1_CS1(1001)+EV1_CS2(1001)+EV1_CS3(1001)+EV1_CS4(1001)),RT_Ta_EV2(1+EV2_CS1(1001)+EV2_CS2(1001)+EV2_CS3(1001):EV2_CS1(1001)+EV2_CS2(1001)+EV2_CS3(1001)+EV2_CS4(1001)),RT_Ta_EV3(1+EV3_CS1(1001)+EV3_CS2(1001)+EV3_CS3(1001):EV3_CS1(1001)+EV3_CS2(1001)+EV3_CS3(1001)+EV3_CS4(1001))];
data_CS1(1001).Tl=[RT_Tl_EV1(1:EV1_CS1(1001)),RT_Tl_EV2(1:EV2_CS1(1001)),RT_Tl_EV3(1:EV3_CS1(1001))];
data_CS2(1001).Tl=[RT_Tl_EV1(1+EV1_CS1(1001):EV1_CS1(1001)+EV1_CS2(1001)),RT_Tl_EV2(1+EV2_CS1(1001):EV2_CS1(1001)+EV2_CS2(1001)),RT_Tl_EV3(1+EV3_CS1(1001):EV3_CS1(1001)+EV3_CS2(1001))];
data_CS3(1001).Tl=[RT_Tl_EV1(1+EV1_CS1(1001)+EV1_CS2(1001):EV1_CS1(1001)+EV1_CS2(1001)+EV1_CS3(1001)),RT_Tl_EV2(1+EV2_CS1(1001)+EV2_CS2(1001):EV2_CS1(1001)+EV2_CS2(1001)+EV2_CS3(1001)),RT_Tl_EV3(1+EV3_CS1(1001)+EV3_CS2(1001):EV3_CS1(1001)+EV3_CS2(1001)+EV3_CS3(1001))];
data_CS4(1001).Tl=[RT_Tl_EV1(1+EV1_CS1(1001)+EV1_CS2(1001)+EV1_CS3(1001):EV1_CS1(1001)+EV1_CS2(1001)+EV1_CS3(1001)+EV1_CS4(1001)),RT_Tl_EV2(1+EV2_CS1(1001)+EV2_CS2(1001)+EV2_CS3(1001):EV2_CS1(1001)+EV2_CS2(1001)+EV2_CS3(1001)+EV2_CS4(1001)),RT_Tl_EV3(1+EV3_CS1(1001)+EV3_CS2(1001)+EV3_CS3(1001):EV3_CS1(1001)+EV3_CS2(1001)+EV3_CS3(1001)+EV3_CS4(1001))];
data_CS1(1001).S0=[RT_S0_EV1(1:EV1_CS1(1001)),RT_S0_EV2(1:EV2_CS1(1001)),RT_S0_EV3(1:EV3_CS1(1001))];
data_CS2(1001).S0=[RT_S0_EV1(1+EV1_CS1(1001):EV1_CS1(1001)+EV1_CS2(1001)),RT_S0_EV2(1+EV2_CS1(1001):EV2_CS1(1001)+EV2_CS2(1001)),RT_S0_EV3(1+EV3_CS1(1001):EV3_CS1(1001)+EV3_CS2(1001))];
data_CS3(1001).S0=[RT_S0_EV1(1+EV1_CS1(1001)+EV1_CS2(1001):EV1_CS1(1001)+EV1_CS2(1001)+EV1_CS3(1001)),RT_S0_EV2(1+EV2_CS1(1001)+EV2_CS2(1001):EV2_CS1(1001)+EV2_CS2(1001)+EV2_CS3(1001)),RT_S0_EV3(1+EV3_CS1(1001)+EV3_CS2(1001):EV3_CS1(1001)+EV3_CS2(1001)+EV3_CS3(1001))];
data_CS4(1001).S0=[RT_S0_EV1(1+EV1_CS1(1001)+EV1_CS2(1001)+EV1_CS3(1001):EV1_CS1(1001)+EV1_CS2(1001)+EV1_CS3(1001)+EV1_CS4(1001)),RT_S0_EV2(1+EV2_CS1(1001)+EV2_CS2(1001)+EV2_CS3(1001):EV2_CS1(1001)+EV2_CS2(1001)+EV2_CS3(1001)+EV2_CS4(1001)),RT_S0_EV3(1+EV3_CS1(1001)+EV3_CS2(1001)+EV3_CS3(1001):EV3_CS1(1001)+EV3_CS2(1001)+EV3_CS3(1001)+EV3_CS4(1001))];
%%%��ʷ�ɵ���Ǳ������%%%
for i=1:1000
    %���վ1
    data_CS1(i).X=zeros(length(data_CS1(i).S0),24);%����ͣ��״̬����
    for j=1:length(data_CS1(i).S0)
        data_CS1(i).X(j,data_CS1(i).Ta(j):data_CS1(i).Tl(j))=1;
    end
    data_CS1(i).Pch=6.6*sum(data_CS1(i).X);%���崢�ܵĳ�繦��
    data_CS1(i).Pdis=6.6*sum(data_CS1(i).X);%���崢�ܵķŵ繦��
    data_CS1(i).Smin(1:23)=32*0.15*sum(data_CS1(i).X(:,1:23))+(32*0.9-32*0.15)*sum(data_CS1(i).X(:,1:23).*(data_CS1(i).X(:,1:23)-data_CS1(i).X(:,2:24)));%���崢�ܵ���С����
    data_CS1(i).Smin(24)=32*0.15*sum(data_CS1(i).X(:,24))+(32*0.9-32*0.15)*sum(data_CS1(i).X(:,24));
    data_CS1(i).Smax=32*0.9*sum(data_CS1(i).X);%���崢�ܵ��������
    data_CS1(i).dS=zeros(1,24);%���崢�ܵ������仯��
    data_CS1(i).dS(1,1)=data_CS1(i).S0*data_CS1(i).X(:,1);
    data_CS1(i).dS(1,2:24)=data_CS1(i).S0*(data_CS1(i).X(:,2:24).*(data_CS1(i).X(:,2:24)-data_CS1(i).X(:,1:23)))-32*0.9*sum(data_CS1(i).X(:,1:23).*(data_CS1(i).X(:,1:23)-data_CS1(i).X(:,2:24)));
    %���վ2
    data_CS2(i).X=zeros(length(data_CS2(i).S0),24);%����ͣ��״̬����
    for j=1:length(data_CS2(i).S0)
        data_CS2(i).X(j,data_CS2(i).Ta(j):data_CS2(i).Tl(j))=1;
    end
    data_CS2(i).Pch=6.6*sum(data_CS2(i).X);%���崢�ܵĳ�繦��
    data_CS2(i).Pdis=6.6*sum(data_CS2(i).X);%���崢�ܵķŵ繦��
    data_CS2(i).Smin(1:23)=32*0.15*sum(data_CS2(i).X(:,1:23))+(32*0.9-32*0.15)*sum(data_CS2(i).X(:,1:23).*(data_CS2(i).X(:,1:23)-data_CS2(i).X(:,2:24)));%���崢�ܵ���С����
    data_CS2(i).Smin(24)=32*0.15*sum(data_CS2(i).X(:,24))+(32*0.9-32*0.15)*sum(data_CS2(i).X(:,24));
    data_CS2(i).Smax=32*0.9*sum(data_CS2(i).X);%���崢�ܵ��������
    data_CS2(i).dS=zeros(1,24);%���崢�ܵ������仯��
    data_CS2(i).dS(1,1)=data_CS2(i).S0*data_CS2(i).X(:,1);
    data_CS2(i).dS(1,2:24)=data_CS2(i).S0*(data_CS2(i).X(:,2:24).*(data_CS2(i).X(:,2:24)-data_CS2(i).X(:,1:23)))-32*0.9*sum(data_CS2(i).X(:,1:23).*(data_CS2(i).X(:,1:23)-data_CS2(i).X(:,2:24)));
    %���վ3
    data_CS3(i).X=zeros(length(data_CS3(i).S0),24);%����ͣ��״̬����
    for j=1:length(data_CS3(i).S0)
        data_CS3(i).X(j,data_CS3(i).Ta(j):data_CS3(i).Tl(j))=1;
    end
    data_CS3(i).Pch=6.6*sum(data_CS3(i).X);%���崢�ܵĳ�繦��
    data_CS3(i).Pdis=6.6*sum(data_CS3(i).X);%���崢�ܵķŵ繦��
    data_CS3(i).Smin(1:23)=32*0.15*sum(data_CS3(i).X(:,1:23))+(32*0.9-32*0.15)*sum(data_CS3(i).X(:,1:23).*(data_CS3(i).X(:,1:23)-data_CS3(i).X(:,2:24)));%���崢�ܵ���С����
    data_CS3(i).Smin(24)=32*0.15*sum(data_CS3(i).X(:,24))+(32*0.9-32*0.15)*sum(data_CS3(i).X(:,24));
    data_CS3(i).Smax=32*0.9*sum(data_CS3(i).X);%���崢�ܵ��������
    data_CS3(i).dS=zeros(1,24);%���崢�ܵ������仯��
    data_CS3(i).dS(1,1)=data_CS3(i).S0*data_CS3(i).X(:,1);
    data_CS3(i).dS(1,2:24)=data_CS3(i).S0*(data_CS3(i).X(:,2:24).*(data_CS3(i).X(:,2:24)-data_CS3(i).X(:,1:23)))-32*0.9*sum(data_CS3(i).X(:,1:23).*(data_CS3(i).X(:,1:23)-data_CS3(i).X(:,2:24)));
    %���վ4
    data_CS4(i).X=zeros(length(data_CS4(i).S0),24);%����ͣ��״̬����
    for j=1:length(data_CS4(i).S0)
        data_CS4(i).X(j,data_CS4(i).Ta(j):data_CS4(i).Tl(j))=1;
    end
    data_CS4(i).Pch=6.6*sum(data_CS4(i).X);%���崢�ܵĳ�繦��
    data_CS4(i).Pdis=6.6*sum(data_CS4(i).X);%���崢�ܵķŵ繦��
    data_CS4(i).Smin(1:23)=32*0.15*sum(data_CS4(i).X(:,1:23))+(32*0.9-32*0.15)*sum(data_CS4(i).X(:,1:23).*(data_CS4(i).X(:,1:23)-data_CS4(i).X(:,2:24)));%���崢�ܵ���С����
    data_CS4(i).Smin(24)=32*0.15*sum(data_CS4(i).X(:,24))+(32*0.9-32*0.15)*sum(data_CS4(i).X(:,24));
    data_CS4(i).Smax=32*0.9*sum(data_CS4(i).X);%���崢�ܵ��������
    data_CS4(i).dS=zeros(1,24);%���崢�ܵ������仯��
    data_CS4(i).dS(1,1)=data_CS4(i).S0*data_CS4(i).X(:,1);
    data_CS4(i).dS(1,2:24)=data_CS4(i).S0*(data_CS4(i).X(:,2:24).*(data_CS4(i).X(:,2:24)-data_CS4(i).X(:,1:23)))-32*0.9*sum(data_CS4(i).X(:,1:23).*(data_CS4(i).X(:,1:23)-data_CS4(i).X(:,2:24)));
end
%%%��ǰ�ɵ���Ǳ��Ԥ��%%%
%���վ1
Pch_CS1=zeros(1000,24);Pdis_CS1=zeros(1000,24);Smin_CS1=zeros(1000,24);Smax_CS1=zeros(1000,24);dS_CS1=zeros(1000,24);
for i=1:1000
    Pch_CS1(i,:)=data_CS1(i).Pch;Pdis_CS1(i,:)=data_CS1(i).Pdis;Smin_CS1(i,:)=data_CS1(i).Smin;Smax_CS1(i,:)=data_CS1(i).Smax;dS_CS1(i,:)=data_CS1(i).dS;
end
Forecast_CS1=[mean(Pch_CS1);mean(Pdis_CS1);mean(Smin_CS1);mean(Smax_CS1);mean(dS_CS1)];
%���վ2
Pch_CS2=zeros(1000,24);Pdis_CS2=zeros(1000,24);Smin_CS2=zeros(1000,24);Smax_CS2=zeros(1000,24);dS_CS2=zeros(1000,24);
for i=1:1000
    Pch_CS2(i,:)=data_CS2(i).Pch;Pdis_CS2(i,:)=data_CS2(i).Pdis;Smin_CS2(i,:)=data_CS2(i).Smin;Smax_CS2(i,:)=data_CS2(i).Smax;dS_CS2(i,:)=data_CS2(i).dS;
end
Forecast_CS2=[mean(Pch_CS2);mean(Pdis_CS2);mean(Smin_CS2);mean(Smax_CS2);mean(dS_CS2)];
%���վ3
Pch_CS3=zeros(1000,24);Pdis_CS3=zeros(1000,24);Smin_CS3=zeros(1000,24);Smax_CS3=zeros(1000,24);dS_CS3=zeros(1000,24);
for i=1:1000
    Pch_CS3(i,:)=data_CS3(i).Pch;Pdis_CS3(i,:)=data_CS3(i).Pdis;Smin_CS3(i,:)=data_CS3(i).Smin;Smax_CS3(i,:)=data_CS3(i).Smax;dS_CS3(i,:)=data_CS3(i).dS;
end
Forecast_CS3=[mean(Pch_CS3);mean(Pdis_CS3);mean(Smin_CS3);mean(Smax_CS3);mean(dS_CS3)];
%���վ4
Pch_CS4=zeros(1000,24);Pdis_CS4=zeros(1000,24);Smin_CS4=zeros(1000,24);Smax_CS4=zeros(1000,24);dS_CS4=zeros(1000,24);
for i=1:1000
    Pch_CS4(i,:)=data_CS4(i).Pch;Pdis_CS4(i,:)=data_CS4(i).Pdis;Smin_CS4(i,:)=data_CS4(i).Smin;Smax_CS4(i,:)=data_CS4(i).Smax;dS_CS4(i,:)=data_CS4(i).dS;
end
Forecast_CS4=[mean(Pch_CS4);mean(Pdis_CS4);mean(Smin_CS4);mean(Smax_CS4);mean(dS_CS4)];
%%%��ǰ�����繦��(���ȳ��ԭ��)%%%
%���վ1
Pch_CS1_disorder=sdpvar(1,24);%���
Pdis_CS1_disorder=zeros(1,24);%�ŵ�
S_CS1_disorder=sdpvar(1,24);%SOC
f_CS1_disorder=Pch_CS1_disorder*[1:24]';%���ȳ��ԭ��Ŀ�꺯��
C_CS1_disorder=[0<=Pch_CS1_disorder<=Forecast_CS1(1,:),
    Forecast_CS1(3,:)<=S_CS1_disorder<=Forecast_CS1(4,:),
    S_CS1_disorder(1)==0.95*Pch_CS1_disorder(1)+Forecast_CS1(5,1),
    S_CS1_disorder(2:24)==S_CS1_disorder(1:23)+0.95*Pch_CS1_disorder(2:24)+Forecast_CS1(5,2:24)];%Լ������
solvesdp(C_CS1_disorder,f_CS1_disorder);
Pch_CS1_disorder=double(Pch_CS1_disorder);
S_CS1_disorder=double(S_CS1_disorder);
%���վ2
Pch_CS2_disorder=sdpvar(1,24);%���
Pdis_CS2_disorder=zeros(1,24);%�ŵ�
S_CS2_disorder=sdpvar(1,24);%SOC
f_CS2_disorder=Pch_CS2_disorder*[1:24]';%���ȳ��ԭ��Ŀ�꺯��
C_CS2_disorder=[0<=Pch_CS2_disorder<=Forecast_CS2(1,:),
    Forecast_CS2(3,:)<=S_CS2_disorder<=Forecast_CS2(4,:),
    S_CS2_disorder(1)==0.95*Pch_CS2_disorder(1)+Forecast_CS2(5,1),
    S_CS2_disorder(2:24)==S_CS2_disorder(1:23)+0.95*Pch_CS2_disorder(2:24)+Forecast_CS2(5,2:24)];%Լ������
solvesdp(C_CS2_disorder,f_CS2_disorder);
Pch_CS2_disorder=double(Pch_CS2_disorder);
S_CS2_disorder=double(S_CS2_disorder);
%���վ3
Pch_CS3_disorder=sdpvar(1,24);%���
Pdis_CS3_disorder=zeros(1,24);%�ŵ�
S_CS3_disorder=sdpvar(1,24);%SOC
f_CS3_disorder=Pch_CS3_disorder*[1:24]';%���ȳ��ԭ��Ŀ�꺯��
C_CS3_disorder=[0<=Pch_CS3_disorder<=Forecast_CS3(1,:),
    Forecast_CS3(3,:)<=S_CS3_disorder<=Forecast_CS3(4,:),
    S_CS3_disorder(1)==0.95*Pch_CS3_disorder(1)+Forecast_CS3(5,1),
    S_CS3_disorder(2:24)==S_CS3_disorder(1:23)+0.95*Pch_CS3_disorder(2:24)+Forecast_CS3(5,2:24)];%Լ������
solvesdp(C_CS3_disorder,f_CS3_disorder);
Pch_CS3_disorder=double(Pch_CS3_disorder);
S_CS3_disorder=double(S_CS3_disorder);
%���վ4
Pch_CS4_disorder=sdpvar(1,24);%���
Pdis_CS4_disorder=zeros(1,24);%�ŵ�
S_CS4_disorder=sdpvar(1,24);%SOC
f_CS4_disorder=Pch_CS4_disorder*[1:24]';%���ȳ��ԭ��Ŀ�꺯��
C_CS4_disorder=[0<=Pch_CS4_disorder<=Forecast_CS4(1,:),
    Forecast_CS4(3,:)<=S_CS4_disorder<=Forecast_CS4(4,:),
    S_CS4_disorder(1)==0.95*Pch_CS4_disorder(1)+Forecast_CS4(5,1),
    S_CS4_disorder(2:24)==S_CS4_disorder(1:23)+0.95*Pch_CS4_disorder(2:24)+Forecast_CS4(5,2:24)];%Լ������
solvesdp(C_CS4_disorder,f_CS4_disorder);
Pch_CS4_disorder=double(Pch_CS4_disorder);
S_CS4_disorder=double(S_CS4_disorder);
%��ͼ���Գ��վ3Ϊ��
figure(1);
hold on
plot(Forecast_CS3(1,:),'b')%��繦�ʱ߽�
plot(-Forecast_CS3(2,:),'g')%�ŵ繦�ʱ߽�
plot(Pch_CS3_disorder,'r.-')
legend('��繦���Ͻ�','�ŵ繦���Ͻ�','ʵ�ʳ�ŵ繦��')
xlabel ʱ��
ylabel ����(kW)
figure(2);
hold on
plot(Forecast_CS3(4,:),'g')%SOC�Ͻ�
plot(Forecast_CS3(3,:),'b')%SOC�½�
plot(S_CS3_disorder,'r.-')
legend('SOC�Ͻ�','SOC�½�','ʵ��SOC')
xlabel ʱ��
ylabel ����(kWh)
%%%ʵʱ�ɵ���Ǳ������(���ýṹ�屣��)%%%
%���վ1
RT_CS1.EV=[data_CS1(1001).Ta;data_CS1(1001).Tl;data_CS1(1001).S0]';
RT_CS1.EV=sortrows(RT_CS1.EV,1);%���ս�վʱ���Ⱥ�˳������
for t=1:96%��ʱ�����ɵ���Ǳ��
    temp1=RT_CS1(1).EV(find(RT_CS1(1).EV(:,1)<=t),:);
    RT_CS1(t).EVset=[temp1];%�õ��綯��������
    [temp2,temp3]=size(RT_CS1(t).EVset);%���󳤶ȣ�Ŀ���ǵõ��綯������������temp2������
    RT_CS1(t).X=zeros(temp2,96);%ͣ��״̬����
    for j=1:temp2
        RT_CS1(t).X(j,RT_CS1(t).EVset(j,1):RT_CS1(t).EVset(j,2))=1;
    end
    RT_CS1(t).EVset=RT_CS1(t).EVset';
    if temp2~=1%ֻ��һ���������Ҫ�������ۣ�����ʹ��sum
        RT_CS1(t).Pch=6.6*sum(RT_CS1(t).X);
        RT_CS1(t).Pdis=6.6*sum(RT_CS1(t).X);
        RT_CS1(t).Smin(1:95)=32*0.15*sum(RT_CS1(t).X(:,1:95))+(32*0.9-32*0.15)*sum(RT_CS1(t).X(:,1:95).*(RT_CS1(t).X(:,1:95)-RT_CS1(t).X(:,2:96)));%���崢�ܵ���С����
        RT_CS1(t).Smin(96)=32*0.9*sum(RT_CS1(t).X(:,96));
        RT_CS1(t).Smax=32*0.9*sum(RT_CS1(t).X);%���崢�ܵ��������
        RT_CS1(t).dS=zeros(1,96);%���崢�ܵ������仯��
        RT_CS1(t).dS(1)=RT_CS1(t).EVset(3,:)*RT_CS1(t).X(:,1);
        RT_CS1(t).dS(2:96)=RT_CS1(t).EVset(3,:)*(RT_CS1(t).X(:,2:96).*(RT_CS1(t).X(:,2:96)-RT_CS1(t).X(:,1:95)))-32*0.9*sum(RT_CS1(t).X(:,1:95).*(RT_CS1(t).X(:,1:95)-RT_CS1(t).X(:,2:96)));
    else
        RT_CS1(t).Pch=6.6*RT_CS1(t).X;
        RT_CS1(t).Pdis=6.6*RT_CS1(t).X;
        RT_CS1(t).Smin(1:95)=32*0.15*RT_CS1(t).X(:,1:95)+(32*0.9-32*0.15)*(RT_CS1(t).X(:,1:95).*(RT_CS1(t).X(:,1:95)-RT_CS1(t).X(:,2:96)));%���崢�ܵ���С����
        RT_CS1(t).Smin(96)=32*0.9*(RT_CS1(t).X(:,96));
        RT_CS1(t).Smax=32*0.9*(RT_CS1(t).X);%���崢�ܵ��������
        RT_CS1(t).dS=zeros(1,96);%���崢�ܵ������仯��
        RT_CS1(t).dS(1)=RT_CS1(t).EVset(3,:)*RT_CS1(t).X(:,1);
        RT_CS1(t).dS(2:96)=RT_CS1(t).EVset(3,:)*(RT_CS1(t).X(:,2:96).*(RT_CS1(t).X(:,2:96)-RT_CS1(t).X(:,1:95)))-32*0.9*(RT_CS1(t).X(:,1:95).*(RT_CS1(t).X(:,1:95)-RT_CS1(t).X(:,2:96)));
    end
end
%���վ2
RT_CS2.EV=[data_CS2(1001).Ta;data_CS2(1001).Tl;data_CS2(1001).S0]';
RT_CS2.EV=sortrows(RT_CS2.EV,1);%���ս�վʱ���Ⱥ�˳������
for t=1:96%��ʱ�����ɵ���Ǳ��
    temp1=RT_CS2(1).EV(find(RT_CS2(1).EV(:,1)<=t),:);
    RT_CS2(t).EVset=[temp1];%�õ��綯��������
    [temp2,temp3]=size(RT_CS2(t).EVset);%���󳤶ȣ�Ŀ���ǵõ��綯������������temp2������
    RT_CS2(t).X=zeros(temp2,96);%ͣ��״̬����
    for j=1:temp2
        RT_CS2(t).X(j,RT_CS2(t).EVset(j,1):RT_CS2(t).EVset(j,2))=1;
    end
    RT_CS2(t).EVset=RT_CS2(t).EVset';
    if temp2~=1%ֻ��һ���������Ҫ�������ۣ�����ʹ��sum
        RT_CS2(t).Pch=6.6*sum(RT_CS2(t).X);
        RT_CS2(t).Pdis=6.6*sum(RT_CS2(t).X);
        RT_CS2(t).Smin(1:95)=32*0.15*sum(RT_CS2(t).X(:,1:95))+(32*0.9-32*0.15)*sum(RT_CS2(t).X(:,1:95).*(RT_CS2(t).X(:,1:95)-RT_CS2(t).X(:,2:96)));%���崢�ܵ���С����
        RT_CS2(t).Smin(96)=32*0.9*sum(RT_CS2(t).X(:,96));
        RT_CS2(t).Smax=32*0.9*sum(RT_CS2(t).X);%���崢�ܵ��������
        RT_CS2(t).dS=zeros(1,96);%���崢�ܵ������仯��
        RT_CS2(t).dS(1)=RT_CS2(t).EVset(3,:)*RT_CS2(t).X(:,1);
        RT_CS2(t).dS(2:96)=RT_CS2(t).EVset(3,:)*(RT_CS2(t).X(:,2:96).*(RT_CS2(t).X(:,2:96)-RT_CS2(t).X(:,1:95)))-32*0.9*sum(RT_CS2(t).X(:,1:95).*(RT_CS2(t).X(:,1:95)-RT_CS2(t).X(:,2:96)));
    else
        RT_CS2(t).Pch=6.6*RT_CS2(t).X;
        RT_CS2(t).Pdis=6.6*RT_CS2(t).X;
        RT_CS2(t).Smin(1:95)=32*0.15*RT_CS2(t).X(:,1:95)+(32*0.9-32*0.15)*(RT_CS2(t).X(:,1:95).*(RT_CS2(t).X(:,1:95)-RT_CS2(t).X(:,2:96)));%���崢�ܵ���С����
        RT_CS2(t).Smin(96)=32*0.9*(RT_CS2(t).X(:,96));
        RT_CS2(t).Smax=32*0.9*(RT_CS2(t).X);%���崢�ܵ��������
        RT_CS2(t).dS=zeros(1,96);%���崢�ܵ������仯��
        RT_CS2(t).dS(1)=RT_CS2(t).EVset(3,:)*RT_CS2(t).X(:,1);
        RT_CS2(t).dS(2:96)=RT_CS2(t).EVset(3,:)*(RT_CS2(t).X(:,2:96).*(RT_CS2(t).X(:,2:96)-RT_CS2(t).X(:,1:95)))-32*0.9*(RT_CS2(t).X(:,1:95).*(RT_CS2(t).X(:,1:95)-RT_CS2(t).X(:,2:96)));
    end
end
%���վ3
RT_CS3.EV=[data_CS3(1001).Ta;data_CS3(1001).Tl;data_CS3(1001).S0]';
RT_CS3.EV=sortrows(RT_CS3.EV,1);%���ս�վʱ���Ⱥ�˳������
for t=1:96%��ʱ�����ɵ���Ǳ��
    temp1=RT_CS3(1).EV(find(RT_CS3(1).EV(:,1)<=t),:);
    RT_CS3(t).EVset=[temp1];%�õ��綯��������
    [temp2,temp3]=size(RT_CS3(t).EVset);%���󳤶ȣ�Ŀ���ǵõ��綯������������temp2������
    RT_CS3(t).X=zeros(temp2,96);%ͣ��״̬����
    for j=1:temp2
        RT_CS3(t).X(j,RT_CS3(t).EVset(j,1):RT_CS3(t).EVset(j,2))=1;
    end
    RT_CS3(t).EVset=RT_CS3(t).EVset';
    if temp2~=1%ֻ��һ���������Ҫ�������ۣ�����ʹ��sum
        RT_CS3(t).Pch=6.6*sum(RT_CS3(t).X);
        RT_CS3(t).Pdis=6.6*sum(RT_CS3(t).X);
        RT_CS3(t).Smin(1:95)=32*0.15*sum(RT_CS3(t).X(:,1:95))+(32*0.9-32*0.15)*sum(RT_CS3(t).X(:,1:95).*(RT_CS3(t).X(:,1:95)-RT_CS3(t).X(:,2:96)));%���崢�ܵ���С����
        RT_CS3(t).Smin(96)=32*0.9*sum(RT_CS3(t).X(:,96));
        RT_CS3(t).Smax=32*0.9*sum(RT_CS3(t).X);%���崢�ܵ��������
        RT_CS3(t).dS=zeros(1,96);%���崢�ܵ������仯��
        RT_CS3(t).dS(1)=RT_CS3(t).EVset(3,:)*RT_CS3(t).X(:,1);
        RT_CS3(t).dS(2:96)=RT_CS3(t).EVset(3,:)*(RT_CS3(t).X(:,2:96).*(RT_CS3(t).X(:,2:96)-RT_CS3(t).X(:,1:95)))-32*0.9*sum(RT_CS3(t).X(:,1:95).*(RT_CS3(t).X(:,1:95)-RT_CS3(t).X(:,2:96)));
    else
        RT_CS3(t).Pch=6.6*RT_CS3(t).X;
        RT_CS3(t).Pdis=6.6*RT_CS3(t).X;
        RT_CS3(t).Smin(1:95)=32*0.15*RT_CS3(t).X(:,1:95)+(32*0.9-32*0.15)*(RT_CS3(t).X(:,1:95).*(RT_CS3(t).X(:,1:95)-RT_CS3(t).X(:,2:96)));%���崢�ܵ���С����
        RT_CS3(t).Smin(96)=32*0.9*(RT_CS3(t).X(:,96));
        RT_CS3(t).Smax=32*0.9*(RT_CS3(t).X);%���崢�ܵ��������
        RT_CS3(t).dS=zeros(1,96);%���崢�ܵ������仯��
        RT_CS3(t).dS(1)=RT_CS3(t).EVset(3,:)*RT_CS3(t).X(:,1);
        RT_CS3(t).dS(2:96)=RT_CS3(t).EVset(3,:)*(RT_CS3(t).X(:,2:96).*(RT_CS3(t).X(:,2:96)-RT_CS3(t).X(:,1:95)))-32*0.9*(RT_CS3(t).X(:,1:95).*(RT_CS3(t).X(:,1:95)-RT_CS3(t).X(:,2:96)));
    end
end
%���վ4
RT_CS4.EV=[data_CS4(1001).Ta;data_CS4(1001).Tl;data_CS4(1001).S0]';
RT_CS4.EV=sortrows(RT_CS4.EV,1);%���ս�վʱ���Ⱥ�˳������
for t=1:96%��ʱ�����ɵ���Ǳ��
    temp1=RT_CS4(1).EV(find(RT_CS4(1).EV(:,1)<=t),:);
    RT_CS4(t).EVset=[temp1];%�õ��綯��������
    [temp2,temp3]=size(RT_CS4(t).EVset);%���󳤶ȣ�Ŀ���ǵõ��綯������������temp2������
    RT_CS4(t).X=zeros(temp2,96);%ͣ��״̬����
    for j=1:temp2
        RT_CS4(t).X(j,RT_CS4(t).EVset(j,1):RT_CS4(t).EVset(j,2))=1;
    end
    RT_CS4(t).EVset=RT_CS4(t).EVset';
    if temp2~=1%ֻ��һ���������Ҫ�������ۣ�����ʹ��sum
        RT_CS4(t).Pch=6.6*sum(RT_CS4(t).X);
        RT_CS4(t).Pdis=6.6*sum(RT_CS4(t).X);
        RT_CS4(t).Smin(1:95)=32*0.15*sum(RT_CS4(t).X(:,1:95))+(32*0.9-32*0.15)*sum(RT_CS4(t).X(:,1:95).*(RT_CS4(t).X(:,1:95)-RT_CS4(t).X(:,2:96)));%���崢�ܵ���С����
        RT_CS4(t).Smin(96)=32*0.9*sum(RT_CS4(t).X(:,96));
        RT_CS4(t).Smax=32*0.9*sum(RT_CS4(t).X);%���崢�ܵ��������
        RT_CS4(t).dS=zeros(1,96);%���崢�ܵ������仯��
        RT_CS4(t).dS(1)=RT_CS4(t).EVset(3,:)*RT_CS4(t).X(:,1);
        RT_CS4(t).dS(2:96)=RT_CS4(t).EVset(3,:)*(RT_CS4(t).X(:,2:96).*(RT_CS4(t).X(:,2:96)-RT_CS4(t).X(:,1:95)))-32*0.9*sum(RT_CS4(t).X(:,1:95).*(RT_CS4(t).X(:,1:95)-RT_CS4(t).X(:,2:96)));
    else
        RT_CS4(t).Pch=6.6*RT_CS4(t).X;
        RT_CS4(t).Pdis=6.6*RT_CS4(t).X;
        RT_CS4(t).Smin(1:95)=32*0.15*RT_CS4(t).X(:,1:95)+(32*0.9-32*0.15)*(RT_CS4(t).X(:,1:95).*(RT_CS4(t).X(:,1:95)-RT_CS4(t).X(:,2:96)));%���崢�ܵ���С����
        RT_CS4(t).Smin(96)=32*0.9*(RT_CS4(t).X(:,96));
        RT_CS4(t).Smax=32*0.9*(RT_CS4(t).X);%���崢�ܵ��������
        RT_CS4(t).dS=zeros(1,96);%���崢�ܵ������仯��
        RT_CS4(t).dS(1)=RT_CS4(t).EVset(3,:)*RT_CS4(t).X(:,1);
        RT_CS4(t).dS(2:96)=RT_CS4(t).EVset(3,:)*(RT_CS4(t).X(:,2:96).*(RT_CS4(t).X(:,2:96)-RT_CS4(t).X(:,1:95)))-32*0.9*(RT_CS4(t).X(:,1:95).*(RT_CS4(t).X(:,1:95)-RT_CS4(t).X(:,2:96)));
    end
end
%��ͼ���Գ��վ3Ϊ��
figure(3)%��24��15����
hold on
plot(RT_CS3(24).Pch,'b')%��繦�ʱ߽�
plot(-RT_CS3(24).Pdis,'g')%�ŵ繦�ʱ߽�
legend('��繦���Ͻ�','�ŵ繦���Ͻ�')
xlabel ʱ��
ylabel ����(kW)
figure(4);
hold on
plot(RT_CS3(24).Smin,'b')%SOC�Ͻ�
plot(RT_CS3(24).Smax,'g')%SOC�½�
legend('SOC�Ͻ�','SOC�½�')
xlabel ʱ��
ylabel ����(kWh)
figure(5)%����ǰ�Ա�
tol=(24-1)/(96-1);%ʱ����
step=1:tol:24;
hold on
plot(Forecast_CS3(1,:),'b')%��繦�ʱ߽�
plot(-Forecast_CS3(2,:),'g')%�ŵ繦�ʱ߽�
plot(step,RT_CS3(96).Pch,'b.-')%��繦�ʱ߽�
plot(step,-RT_CS3(96).Pdis,'g.-')%�ŵ繦�ʱ߽�
legend('��ǰ��繦���Ͻ�','��ǰ�ŵ繦���Ͻ�')
xlabel ʱ��
ylabel ����(kW)
figure(6);
hold on
plot(Forecast_CS3(4,:),'b')%SOC�Ͻ�
plot(Forecast_CS3(3,:),'g')%SOC�½�
plot(step,RT_CS3(96).Smin,'b.-')%SOC�Ͻ�
plot(step,RT_CS3(96).Smax,'g.-')%SOC�½�
legend('��ǰSOC�Ͻ�','��ǰSOC�½�')
xlabel ʱ��
ylabel ����(kWh)
%%%���ݼ���%%%
save('data_DA_potential','Forecast_CS1','Forecast_CS2','Forecast_CS3','Forecast_CS4');%��ǰ�ɵ���Ǳ��
save('data_RT_potential','RT_CS1','RT_CS2','RT_CS3','RT_CS4');%ʵʱ�ɵ���Ǳ��
save('data_disorder','Pch_CS1_disorder','Pch_CS2_disorder','Pch_CS3_disorder','Pch_CS4_disorder');%����������
save('data_EV','data_CS1','data_CS2','data_CS3','data_CS4');%�綯������������
%% �ڶ��ڣ�������Ͷ��
clear
clc
close all
load data_disorder
Pch=[Pch_CS1_disorder;Pch_CS2_disorder;Pch_CS3_disorder;Pch_CS4_disorder];%���վ��繦��
%�г���������
Loadcurve=[0.955391944564747,0.978345604157644,1,0.995019488956258,0.972932005197055,0.970333477695972,0.930489389346037,0.890428757037679,0.902771762667822,0.941966219142486,0.911000433087917,0.862061498484192,0.840190558683413,0.831095712429623,0.756604590731919,0.671719359029883,0.611520138588133,0.582936336076224,0.572542226071893,0.574707665656128,0.587267215244695,0.644218276310091,0.755521870939801,0.884798614118666];
PL_base=[5.704;5.705;5.631;6.518;4.890;5.705;5.847]*1000;%���ɷֲ�
PL=PL_base*Loadcurve;%��������(�������ߴ�08:00��ʼ���𣬼���9��ʱ��)
Pf=sdpvar(7,24);%���߹���
Pf(1,:)=PL(1,:)+Pch(1,:);Pf(2,:)=PL(2,:);Pf(3,:)=PL(3,:);Pf(4,:)=PL(4,:)+Pch(2,:);Pf(5,:)=PL(5,:)+Pch(3,:);Pf(6,:)=PL(6,:);Pf(7,:)=PL(7,:)+Pch(4,:);%���߹������
Pg=sdpvar(10,24);%�����̷ֶε���
Pg_step=1000*[20,5,3,2,2,2,2,2,2,inf]';%��������
Price_DSO=[3:12]'*0.1;%�ֶε��
Obj=sum(sum((Price_DSO*ones(1,24)).*Pg));%Ŀ��Ϊ�õ������С
Constraint=[0<=Pg<=Pg_step*ones(1,24),sum(Pg)==sum(Pf)];%Լ������
solvesdp(Constraint,Obj);%������Թ滮����
Pg=double(Pg);%���������
Pf=double(Pf);%���߹���
isPg=(Pg>0);%Ϊ�˼�������ۣ����㷢����ֶ�ѡ�����
DLMP=sum(isPg)/10+0.2;%�����ۼ���
%��ͼ
figure(1)%�ڵ�߼ʵ��
stairs(DLMP);
xlabel ʱ��
ylabel ���(Ԫ/kWh)
ylim([0.3,1.3])
figure(2)%��������
hold on
plot(sum(PL)/1000);
plot(sum(Pf)/1000,'r.-');
xlabel ʱ��
ylabel ����(MW)
legend('��������','�����縺��')
Cost=sum(sum(Pch).*DLMP);%���õ����
result_disorder.Cost=Cost;result_disorder.DLMP=DLMP;result_disorder.Pf=Pf;result_disorder.Pg=Pg;%�������
save('result_disorder','result_disorder');
%% �����ڣ���ǰͶ�����
clear
clc
close all
load data_DA_potential
Pchmax=[Forecast_CS1(1,:);Forecast_CS2(1,:);Forecast_CS3(1,:);Forecast_CS4(1,:)];%��繦������
Pdismax=[Forecast_CS1(2,:);Forecast_CS2(2,:);Forecast_CS3(2,:);Forecast_CS4(2,:)];%�ŵ繦������
Smin=[Forecast_CS1(3,:);Forecast_CS2(3,:);Forecast_CS3(3,:);Forecast_CS4(3,:)];%SOC����
Smax=[Forecast_CS1(4,:);Forecast_CS2(4,:);Forecast_CS3(4,:);Forecast_CS4(4,:)];%SOC����
dS=[Forecast_CS1(5,:);Forecast_CS2(5,:);Forecast_CS3(5,:);Forecast_CS4(5,:)];%SOC�仯
%���վͶ������
Pch=sdpvar(4,24);%���վ��繦��
Pdis=sdpvar(4,24);%���վ�ŵ繦��
S=sdpvar(4,24);%���վ���崢��SOC
Ccs=[0<=Pch<=Pchmax,0<=Pdis<=Pdismax,Smin<=S<=Smax,S(:,1)==0.95*Pch(:,1)-Pdis(:,1)/0.95+dS(:,1),S(:,2:24)==S(:,1:23)+0.95*Pch(:,2:24)-Pdis(:,2:24)/0.95+dS(:,2:24)];%Լ������
% %�г���������
Loadcurve=[0.955391944564747,0.978345604157644,1,0.995019488956258,0.972932005197055,0.970333477695972,0.930489389346037,0.890428757037679,0.902771762667822,0.941966219142486,0.911000433087917,0.862061498484192,0.840190558683413,0.831095712429623,0.756604590731919,0.671719359029883,0.611520138588133,0.582936336076224,0.572542226071893,0.574707665656128,0.587267215244695,0.644218276310091,0.755521870939801,0.884798614118666];
PL_base=[5.704;5.705;5.631;6.518;4.890;5.705;5.847]*1000;%���ɷֲ�
PL=PL_base*Loadcurve;%��������
Pf=sdpvar(7,24);%���߹���
Pf(1,:)=PL(1,:)+Pch(1,:)-Pdis(1,:);Pf(2,:)=PL(2,:);Pf(3,:)=PL(3,:);Pf(4,:)=PL(4,:)+Pch(2,:)-Pdis(2,:);Pf(5,:)=PL(5,:)+Pch(3,:)-Pdis(3,:);Pf(6,:)=PL(6,:);Pf(7,:)=PL(7,:)+Pch(4,:)-Pdis(4,:);%���߹������
Pf_limit=1000*[40,40,40,40,40,40,40]';%���߹�������
Pg=sdpvar(10,24);%�����̷ֶε���
Pg_step=1000*[20,5,3,2,2,2,2,2,2,100]';%��������
Price_DSO=[3:12]'*0.1;%�ֶε��
Lagrant_G_left=sdpvar(10,24);%�����̵����½�
Lagrant_G_right=sdpvar(10,24);%�����̵����Ͻ�
b_Lagrant_G_left=binvar(10,24);%�����̵����½粼������
b_Lagrant_G_right=binvar(10,24);%�����̵����Ͻ粼������
Lagrant_L_left=sdpvar(7,24);%��·�����½�
Lagrant_L_right=sdpvar(7,24);%��·�����Ͻ�
b_Lagrant_L_left=binvar(7,24);%��·�����Ͻ粼������
b_Lagrant_L_right=binvar(7,24);%��·�����½粼������
Lagrant_G=sdpvar(1,24);%ƽ��ڵ���
DLMP=sdpvar(7,24);%DLMP
Ckkt=[sum(Pg)==sum(Pf),
    0<=Pg<=Pg_step*ones(1,24),
    Price_DSO*ones(1,24)-Lagrant_G_left+Lagrant_G_right-ones(10,1)*Lagrant_G==0,
    DLMP==ones(7,1)*Lagrant_G+Lagrant_L_right-Lagrant_L_left,
    Pg<=1E6*b_Lagrant_G_left,0<=Lagrant_G_left<=(1-b_Lagrant_G_left),
    Pg_step*ones(1,24)-Pg<=1E6*b_Lagrant_G_right,0<=Lagrant_G_right<=(1-b_Lagrant_G_right),
    -Pf_limit*ones(1,24)<=Pf<=Pf_limit*ones(1,24),
    Pf+Pf_limit*ones(1,24)<=1E6*b_Lagrant_L_left,0<=Lagrant_L_left<=(1-b_Lagrant_L_left),
    -Pf+Pf_limit*ones(1,24)<=1E6*b_Lagrant_L_right,0<=Lagrant_L_right<=(1-b_Lagrant_L_right)];%�г����������KKT����
Obj=sum(sum((Price_DSO*ones(1,24)).*Pg))-sum(sum(DLMP.*PL))+sum(sum((Pg_step*ones(1,24)).*Lagrant_G_right))+sum(sum((Lagrant_L_left+Lagrant_L_right).*(Pf_limit*ones(1,24))));
%�������
C=[Ccs,Ckkt];%Լ������
ops=sdpsettings('solver','gurobi');
result=solvesdp(C,Obj,ops);
DLMP=double(DLMP);
Lagrant_G=double(Lagrant_G);
Pf=double(Pf);
Pch=double(Pch);
Pdis=double(Pdis);
Pg=double(Pg);
S=double(S);
figure(1)%�ڵ�߼ʵ��
stairs(mean(DLMP));
xlabel ʱ��
ylabel ���(Ԫ/kWh)
ylim([0.3,1.3])
figure(2)%��������
plot(sum(PL),'b--');
hold on
plot(sum(Pg),'r.-')
xlabel ʱ��
ylabel ����(MW)
legend('��������','�����縺��')
figure(3)%��������崢�ܳ�ŵ繦��
plot(sum(Pchmax));
hold on
plot(-sum(Pdismax),'g');
plot(sum(Pch-Pdis),'r.-')
xlabel ʱ��
ylabel ����(kW)
legend('��ǰ��繦���Ͻ�','��ǰ�ŵ繦���Ͻ�','�����ŵ繦��')
figure(4)%��������崢��SOC�仯
plot(sum(Smin),'g');
hold on
plot(sum(Smax));
plot(sum(S),'r.-');
xlabel ʱ��
ylabel ����(kWh)
legend('SOC�Ͻ�','SOC�½�','�����ŵ�SOC')
Cost=sum(sum((Pf-PL).*DLMP));%���õ����
Cost_b=sum(sum(((Pf-PL)>0).*(Pf-PL).*DLMP));%�ܹ������
Cost_s=-sum(sum(((Pf-PL)<0).*(Pf-PL).*DLMP));%���۵�����
result_order.Pdis=Pdis;result_order.Pch=Pch;result_order.Cost_b=Cost_b;result_order.Cost_s=Cost_s;result_order.Cost=Cost;result_order.DLMP=DLMP;result_order.Pf=Pf;result_order.Pg=Pg;%�������
save('result_order','result_order');
%% ���Ľڣ�ʵʱͶ�����
clear
clc
close all
load data_RT_potential
load result_order
Pch1=RT_CS1(96).Pch;Pdis1=RT_CS1(96).Pdis;Smin1=RT_CS1(96).Smin;Smax1=RT_CS1(96).Smax;dS1=RT_CS1(96).dS;%ȡ���һ��ʱ�̵Ŀɵ���Ǳ��
Pch2=RT_CS2(96).Pch;Pdis2=RT_CS2(96).Pdis;Smin2=RT_CS2(96).Smin;Smax2=RT_CS2(96).Smax;dS2=RT_CS2(96).dS;
Pch3=RT_CS3(96).Pch;Pdis3=RT_CS3(96).Pdis;Smin3=RT_CS3(96).Smin;Smax3=RT_CS3(96).Smax;dS3=RT_CS3(96).dS;
Pch4=RT_CS4(96).Pch;Pdis4=RT_CS4(96).Pdis;Smin4=RT_CS4(96).Smin;Smax4=RT_CS4(96).Smax;dS4=RT_CS4(96).dS;
Pchmax=[Pch1;Pch2;Pch3;Pch4];
Pdismax=[Pdis1;Pdis2;Pdis3;Pdis4];
Smin=[Smin1;Smin2;Smin3;Smin4];
Smax=[Smax1;Smax2;Smax3;Smax4];
dS=[dS1;dS2;dS3;dS4];
Link=zeros(24,96);%ʱ�λ������(��ǰ1h����Ϊʵʱ15min)
for i=1:24
    Link(i,4*i-3:4*i)=1;
end
DA_Ps=result_order.Pch-result_order.Pdis;%��ǰ��������
DA_Ps=DA_Ps*Link;%ʱ�λ���
DA_Pf=result_order.Pf*Link;%��ǰ�ܹ���
%���վͶ������
Pch=sdpvar(4,96);%���վ��繦��
Pdis=sdpvar(4,96);%���վ�ŵ繦��
S=sdpvar(4,96);%���վ���崢��SOC
Ccs=[0<=Pch<=Pchmax,0<=Pdis<=Pdismax,Smin<=S<=Smax,S(:,1)==0.95*0.25*Pch(:,1)-0.25*Pdis(:,1)/0.95+dS(:,1),S(:,2:96)==S(:,1:95)+0.25*0.95*Pch(:,2:96)-0.25*Pdis(:,2:96)/0.95+dS(:,2:96)];%Լ������
% %�г���������
Loadcurve=[0.955391944564747,0.978345604157644,1,0.995019488956258,0.972932005197055,0.970333477695972,0.930489389346037,0.890428757037679,0.902771762667822,0.941966219142486,0.911000433087917,0.862061498484192,0.840190558683413,0.831095712429623,0.756604590731919,0.671719359029883,0.611520138588133,0.582936336076224,0.572542226071893,0.574707665656128,0.587267215244695,0.644218276310091,0.755521870939801,0.884798614118666];
Loadcurve=Loadcurve*Link;%����96��ʱ��
PL_base=[5.704;5.705;5.631;6.518;4.890;5.705;5.847]*1000;%���ɷֲ�
PL=PL_base*Loadcurve;%��������
Pf=sdpvar(7,96);%���߹���
Pf(1,:)=PL(1,:)+Pch(1,:)-Pdis(1,:);Pf(2,:)=PL(2,:);Pf(3,:)=PL(3,:);Pf(4,:)=PL(4,:)+Pch(2,:)-Pdis(2,:);Pf(5,:)=PL(5,:)+Pch(3,:)-Pdis(3,:);Pf(6,:)=PL(6,:);Pf(7,:)=PL(7,:)+Pch(4,:)-Pdis(4,:);%���߹������
Pf_limit=1000*[40,40,40,40,40,40,40]';%���߹�������
Pg=sdpvar(10,96);%�����̷ֶε���
Pg_step=1000*[20,5,3,2,2,2,2,2,2,100]';%��������
Price_DSO=[3:12]'*0.1;%�ֶε��
Lagrant_G_left=sdpvar(10,96);%�����̵����½�
Lagrant_G_right=sdpvar(10,96);%�����̵����Ͻ�
b_Lagrant_G_left=binvar(10,96);%�����̵����½粼������
b_Lagrant_G_right=binvar(10,96);%�����̵����Ͻ粼������
Lagrant_L_left=sdpvar(7,96);%��·�����½�
Lagrant_L_right=sdpvar(7,96);%��·�����Ͻ�
b_Lagrant_L_left=binvar(7,96);%��·�����Ͻ粼������
b_Lagrant_L_right=binvar(7,96);%��·�����½粼������
Lagrant_G=sdpvar(1,96);%ƽ��ڵ���
DLMP=sdpvar(7,96);%DLMP
Ckkt=[sum(Pg)==sum(Pf),
    0<=Pg<=Pg_step*ones(1,96),
    0.25*Price_DSO*ones(1,96)-Lagrant_G_left+Lagrant_G_right-ones(10,1)*Lagrant_G==0,
    DLMP==ones(7,1)*Lagrant_G+Lagrant_L_right-Lagrant_L_left,
    Pg<=1E6*b_Lagrant_G_left,0<=Lagrant_G_left<=(1-b_Lagrant_G_left),
    Pg_step*ones(1,96)-Pg<=1E6*b_Lagrant_G_right,0<=Lagrant_G_right<=(1-b_Lagrant_G_right),
    -Pf_limit*ones(1,96)<=Pf<=Pf_limit*ones(1,96),
    Pf+Pf_limit*ones(1,96)<=1E6*b_Lagrant_L_left,0<=Lagrant_L_left<=(1-b_Lagrant_L_left),
    -Pf+Pf_limit*ones(1,96)<=1E6*b_Lagrant_L_right,0<=Lagrant_L_right<=(1-b_Lagrant_L_right)];%�г����������KKT����
Obj=0.25*sum(sum((Price_DSO*ones(1,96)).*Pg))-sum(sum(DLMP.*PL))+sum(sum((Pg_step*ones(1,96)).*Lagrant_G_right))+sum(sum((Lagrant_L_left+Lagrant_L_right).*(Pf_limit*ones(1,96))))+0.1*0.25*sum(sum(abs(Pch-Pdis-DA_Ps)));
%�������
C=[Ccs,Ckkt];%Լ������
result=solvesdp(C,Obj);
DLMP=double(DLMP);
Lagrant_G=double(Lagrant_G);
Pf=double(Pf);
Pch=double(Pch);
Pdis=double(Pdis);
Pg=double(Pg);
S=double(S);
double(Obj)
figure(1)%ʵʱ��ǰ����Ա�
hold on
plot(sum(Pf),'r.-')
plot(sum(DA_Pf),'b--');
Cost=sum(sum((Pf-PL).*DLMP));%���õ����
Cost_b=sum(sum(((Pf-PL)>0).*(Pf-PL).*DLMP))%�ܹ������
Cost_s=-sum(sum(((Pf-PL)<0).*(Pf-PL).*DLMP))%���۵�����
Cost_r=0.1*0.25*sum(sum(abs(Pch-Pdis-DA_Ps)))%ƽ�����
MAPE=100*sum(abs((sum(DA_Pf)-sum(Pf))./sum(DA_Pf))/96);%ƽ���ٷ����
result_order_RT.Pdis=Pdis;result_order_RT.Pch=Pch;result_order_RT.Cost_b=Cost_b;result_order_RT.Cost_s=Cost_s;result_order_RT.Cost=Cost;result_order_RT.DLMP=DLMP;result_order_RT.Pf=Pf;result_order_RT.Pg=Pg;%�������
save('result_order_RT','result_order_RT');
%% ����ڣ�ʵʱ��ŵ繦�ʷ���
clear
clc
close all
load result_order_RT
load data_RT_potential
%���վ1
[m,n]=size(RT_CS1(1).EV);%��������
X_CS1=zeros(m,96);
for i=1:m
    X_CS1(i,RT_CS1(1).EV(i,1):RT_CS1(1).EV(i,2))=1;%ͣ��״̬����
end
pch_CS1=sdpvar(m,96);%���
pdis_CS1=sdpvar(m,96);%�ŵ�
s_CS1=sdpvar(m,96);%SOC
C_CS1=[0<=pch_CS1<=6.6*X_CS1,
    0<=pdis_CS1<=6.6*X_CS1,
    32*0.15<=s_CS1<=32*0.9,
    s_CS1(:,1)==RT_CS1(1).EV(:,3)+0.95*0.25*pch_CS1(:,1)-0.25*pdis_CS1(:,1)/0.95,
    s_CS1(:,2:96)==s_CS1(:,1:95)+0.95*0.25*pch_CS1(:,2:96)-0.25*pdis_CS1(:,2:96)/0.95,
    0.25*sum(0.95*pch_CS1-pdis_CS1/0.95,2)==32*0.9-RT_CS1(1).EV(:,3)];%Լ������
a_CS1=sdpvar(1,96);%���������������ڴ治��
b_CS1=sdpvar(1,96);%���������������ڴ治��
C_CS1=[C_CS1,a_CS1==sum(pch_CS1)-result_order_RT.Pch(1,:)];
C_CS1=[C_CS1,b_CS1==sum(pdis_CS1)-result_order_RT.Pdis(1,:)];
f_CS1=sum(a_CS1.^2+b_CS1.^2);
f_CS1=sum(abs(a_CS1)+abs(b_CS1));
solvesdp(C_CS1,f_CS1);
pch_CS1=double(pch_CS1);
pdis_CS1=double(pdis_CS1);
%���վ2
[m,n]=size(RT_CS2(1).EV);%��������
X_CS2=zeros(m,96);
for i=1:m
    X_CS2(i,RT_CS2(1).EV(i,1):RT_CS2(1).EV(i,2))=1;%ͣ��״̬����
end
pch_CS2=sdpvar(m,96);%���
pdis_CS2=sdpvar(m,96);%�ŵ�
s_CS2=sdpvar(m,96);%SOC
C_CS2=[0<=pch_CS2<=6.6*X_CS2,
    0<=pdis_CS2<=6.6*X_CS2,
    32*0.15<=s_CS2<=32*0.9,
    s_CS2(:,1)==RT_CS2(1).EV(:,3)+0.95*0.25*pch_CS2(:,1)-0.25*pdis_CS2(:,1)/0.95,
    s_CS2(:,2:96)==s_CS2(:,1:95)+0.95*0.25*pch_CS2(:,2:96)-0.25*pdis_CS2(:,2:96)/0.95,
    0.25*sum(0.95*pch_CS2-pdis_CS2/0.95,2)==32*0.9-RT_CS2(1).EV(:,3)];%Լ������
a_CS2=sdpvar(1,96);%���������������ڴ治��
b_CS2=sdpvar(1,96);%���������������ڴ治��
C_CS2=[C_CS2,a_CS2==sum(pch_CS2)-result_order_RT.Pch(2,:)];
C_CS2=[C_CS2,b_CS2==sum(pdis_CS2)-result_order_RT.Pdis(2,:)];
f_CS2=sum(a_CS2.^2+b_CS2.^2);
f_CS2=sum(abs(a_CS2)+abs(b_CS2));
solvesdp(C_CS2,f_CS2);
pch_CS2=double(pch_CS2);
pdis_CS2=double(pdis_CS2);
%���վ3
[m,n]=size(RT_CS3(1).EV);%��������
X_CS3=zeros(m,96);
for i=1:m
    X_CS3(i,RT_CS3(1).EV(i,1):RT_CS3(1).EV(i,2))=1;%ͣ��״̬����
end
pch_CS3=sdpvar(m,96);%���
pdis_CS3=sdpvar(m,96);%�ŵ�
s_CS3=sdpvar(m,96);%SOC
C_CS3=[0<=pch_CS3<=6.6*X_CS3,
    0<=pdis_CS3<=6.6*X_CS3,
    32*0.15<=s_CS3<=32*0.9,
    s_CS3(:,1)==RT_CS3(1).EV(:,3)+0.95*0.25*pch_CS3(:,1)-0.25*pdis_CS3(:,1)/0.95,
    s_CS3(:,2:96)==s_CS3(:,1:95)+0.95*0.25*pch_CS3(:,2:96)-0.25*pdis_CS3(:,2:96)/0.95,
    0.25*sum(0.95*pch_CS3-pdis_CS3/0.95,2)==32*0.9-RT_CS3(1).EV(:,3)];%Լ������
a_CS3=sdpvar(1,96);%���������������ڴ治��
b_CS3=sdpvar(1,96);%���������������ڴ治��
C_CS3=[C_CS3,a_CS3==sum(pch_CS3)-result_order_RT.Pch(3,:)];
C_CS3=[C_CS3,b_CS3==sum(pdis_CS3)-result_order_RT.Pdis(3,:)];
f_CS3=sum(a_CS3.^2+b_CS3.^2);
f_CS3=sum(abs(a_CS3)+abs(b_CS3));
solvesdp(C_CS3,f_CS3);
pch_CS3=double(pch_CS3);
pdis_CS3=double(pdis_CS3);
%���վ4
[m,n]=size(RT_CS4(1).EV);%��������
X_CS4=zeros(m,96);
for i=1:m
    X_CS4(i,RT_CS4(1).EV(i,1):RT_CS4(1).EV(i,2))=1;%ͣ��״̬����
end
pch_CS4=sdpvar(m,96);%���
pdis_CS4=sdpvar(m,96);%�ŵ�
s_CS4=sdpvar(m,96);%SOC
C_CS4=[0<=pch_CS4<=6.6*X_CS4,
    0<=pdis_CS4<=6.6*X_CS4,
    32*0.15<=s_CS4<=32*0.9,
    s_CS4(:,1)==RT_CS4(1).EV(:,3)+0.95*0.25*pch_CS4(:,1)-0.25*pdis_CS4(:,1)/0.95,
    s_CS4(:,2:96)==s_CS4(:,1:95)+0.95*0.25*pch_CS4(:,2:96)-0.25*pdis_CS4(:,2:96)/0.95,
    0.25*sum(0.95*pch_CS4-pdis_CS4/0.95,2)==32*0.9-RT_CS4(1).EV(:,3)];%Լ������
a_CS4=sdpvar(1,96);%���������������ڴ治��
b_CS4=sdpvar(1,96);%���������������ڴ治��
C_CS4=[C_CS4,a_CS4==sum(pch_CS4)-result_order_RT.Pch(4,:)];
C_CS4=[C_CS4,b_CS4==sum(pdis_CS4)-result_order_RT.Pdis(4,:)];
f_CS4=sum(a_CS4.^2+b_CS4.^2);
f_CS4=sum(abs(a_CS4)+abs(b_CS4));
solvesdp(C_CS4,f_CS4);
pch_CS4=double(pch_CS4);
pdis_CS4=double(pdis_CS4);
%% �����ڣ����������ȷ������ŵ粹��ϵ����
clear
clc
close all
%�ŵ粹��ϵ��
load data_DA_potential
Pchmax=[Forecast_CS1(1,:);Forecast_CS2(1,:);Forecast_CS3(1,:);Forecast_CS4(1,:)];%��繦������
Pdismax=[Forecast_CS1(2,:);Forecast_CS2(2,:);Forecast_CS3(2,:);Forecast_CS4(2,:)];%�ŵ繦������
Smin=[Forecast_CS1(3,:);Forecast_CS2(3,:);Forecast_CS3(3,:);Forecast_CS4(3,:)];%SOC����
Smax=[Forecast_CS1(4,:);Forecast_CS2(4,:);Forecast_CS3(4,:);Forecast_CS4(4,:)];%SOC����
dS=[Forecast_CS1(5,:);Forecast_CS2(5,:);Forecast_CS3(5,:);Forecast_CS4(5,:)];%SOC�仯
for i=1:51%���Բ�ͬ�ŵ粹��ϵ��
    co_dis=1+(i-1)*0.1;%�ŵ粹��ϵ��
    data_dis(i)=co_dis;
    %���վͶ������
    Pch=sdpvar(4,24);%���վ��繦��
    Pdis=sdpvar(4,24);%���վ�ŵ繦��
    S=sdpvar(4,24);%���վ���崢��SOC
    Ccs=[0<=Pch<=Pchmax,0<=Pdis<=Pdismax,Smin<=S<=Smax,S(:,1)==0.95*Pch(:,1)-Pdis(:,1)/(0.95/co_dis)+dS(:,1),S(:,2:24)==S(:,1:23)+0.95*Pch(:,2:24)-Pdis(:,2:24)/(0.95/co_dis)+dS(:,2:24)];%Լ������
    % %�г���������
    Loadcurve=[0.955391944564747,0.978345604157644,1,0.995019488956258,0.972932005197055,0.970333477695972,0.930489389346037,0.890428757037679,0.902771762667822,0.941966219142486,0.911000433087917,0.862061498484192,0.840190558683413,0.831095712429623,0.756604590731919,0.671719359029883,0.611520138588133,0.582936336076224,0.572542226071893,0.574707665656128,0.587267215244695,0.644218276310091,0.755521870939801,0.884798614118666];
    PL_base=[5.704;5.705;5.631;6.518;4.890;5.705;5.847]*1000;%���ɷֲ�
    PL=PL_base*Loadcurve;%��������
    Pf=sdpvar(7,24);%���߹���
    Pf(1,:)=PL(1,:)+Pch(1,:)-Pdis(1,:);Pf(2,:)=PL(2,:);Pf(3,:)=PL(3,:);Pf(4,:)=PL(4,:)+Pch(2,:)-Pdis(2,:);Pf(5,:)=PL(5,:)+Pch(3,:)-Pdis(3,:);Pf(6,:)=PL(6,:);Pf(7,:)=PL(7,:)+Pch(4,:)-Pdis(4,:);%���߹������
    Pf_limit=1000*[40,40,40,40,40,40,40]';%���߹�������
    Pg=sdpvar(10,24);%�����̷ֶε���
    Pg_step=1000*[20,5,3,2,2,2,2,2,2,100]';%��������
    Price_DSO=[3:12]'*0.1;%�ֶε��
    Lagrant_G_left=sdpvar(10,24);%�����̵����½�
    Lagrant_G_right=sdpvar(10,24);%�����̵����Ͻ�
    b_Lagrant_G_left=binvar(10,24);%�����̵����½粼������
    b_Lagrant_G_right=binvar(10,24);%�����̵����Ͻ粼������
    Lagrant_L_left=sdpvar(7,24);%��·�����½�
    Lagrant_L_right=sdpvar(7,24);%��·�����Ͻ�
    b_Lagrant_L_left=binvar(7,24);%��·�����Ͻ粼������
    b_Lagrant_L_right=binvar(7,24);%��·�����½粼������
    Lagrant_G=sdpvar(1,24);%ƽ��ڵ���
    DLMP=sdpvar(7,24);%DLMP
    Ckkt=[sum(Pg)==sum(Pf),
        0<=Pg<=Pg_step*ones(1,24),
        Price_DSO*ones(1,24)-Lagrant_G_left+Lagrant_G_right-ones(10,1)*Lagrant_G==0,
        DLMP==ones(7,1)*Lagrant_G+Lagrant_L_right-Lagrant_L_left,
        Pg<=1E6*b_Lagrant_G_left,0<=Lagrant_G_left<=(1-b_Lagrant_G_left),
        Pg_step*ones(1,24)-Pg<=1E6*b_Lagrant_G_right,0<=Lagrant_G_right<=(1-b_Lagrant_G_right),
        -Pf_limit*ones(1,24)<=Pf<=Pf_limit*ones(1,24),
        Pf+Pf_limit*ones(1,24)<=1E6*b_Lagrant_L_left,0<=Lagrant_L_left<=(1-b_Lagrant_L_left),
        -Pf+Pf_limit*ones(1,24)<=1E6*b_Lagrant_L_right,0<=Lagrant_L_right<=(1-b_Lagrant_L_right)];%�г����������KKT����
    Obj=sum(sum((Price_DSO*ones(1,24)).*Pg))-sum(sum(DLMP.*PL))+sum(sum((Pg_step*ones(1,24)).*Lagrant_G_right))+sum(sum((Lagrant_L_left+Lagrant_L_right).*(Pf_limit*ones(1,24))));
    %�������
    C=[Ccs,Ckkt];%Լ������
    ops=sdpsettings('solver','gurobi');
    result=solvesdp(C,Obj,ops);
    DLMP=double(DLMP);
    Lagrant_G=double(Lagrant_G);
    Pf=double(Pf);
    Pch=double(Pch);
    Pdis=double(Pdis);
    Pg=double(Pg);
    S=double(S);
    Cost=sum(sum((Pf-PL).*DLMP));%���õ����
    Cost_b=sum(sum(((Pf-PL)>0).*(Pf-PL).*DLMP));%�ܹ������
    Cost_s=-sum(sum(((Pf-PL)<0).*(Pf-PL).*DLMP));%���۵�����
    VD=max(sum(Pf))-min(sum(Pf));%�������Ȳ�
    data_cost(i)=Cost;
    data_cost_b(i)=Cost_b;
    data_cost_s(i)=Cost_s;
    data_VD(i)=VD;
end
%% ���߽ڣ���ǰ�ɵ���Ǳ��Ԥ���㷨�Ա�
%�㷨1������
clear
clc
close all
load data_EV
X_CS1=data_CS1.X;X_CS2=data_CS2.X;X_CS3=data_CS3.X;X_CS4=data_CS4.X;
S0_CS1=data_CS1.S0;S0_CS2=data_CS2.S0;S0_CS3=data_CS3.S0;S0_CS4=data_CS4.S0;
m1=length(S0_CS1);m2=length(S0_CS2);m3=length(S0_CS3);m4=length(S0_CS4);
%���վͶ������
pch_CS1=sdpvar(m1,24);pdis_CS1=sdpvar(m1,24);%�綯��������ĳ�ŵ繦��
pch_CS2=sdpvar(m2,24);pdis_CS2=sdpvar(m2,24);
pch_CS3=sdpvar(m3,24);pdis_CS3=sdpvar(m3,24);
pch_CS4=sdpvar(m4,24);pdis_CS4=sdpvar(m4,24);
s_CS1=sdpvar(m1,24);s_CS2=sdpvar(m2,24);s_CS3=sdpvar(m3,24);s_CS4=sdpvar(m4,24);%�綯���������SOC
Pch=sdpvar(4,24);%���վ��繦��
Pdis=sdpvar(4,24);%���վ�ŵ繦��
Ccs=[Pch(1,:)==sum(pch_CS1),Pch(2,:)==sum(pch_CS2),Pch(3,:)==sum(pch_CS3),Pch(4,:)==sum(pch_CS4),
    Pdis(1,:)==sum(pdis_CS1),Pdis(2,:)==sum(pdis_CS2),Pdis(3,:)==sum(pdis_CS3),Pdis(4,:)==sum(pdis_CS4),
    0<=pch_CS1<=6.6*X_CS1,0<=pch_CS2<=6.6*X_CS2,0<=pch_CS3<=6.6*X_CS3,0<=pch_CS4<=6.6*X_CS4,
    0<=pdis_CS1<=6.6*X_CS1,0<=pdis_CS2<=6.6*X_CS2,0<=pdis_CS3<=6.6*X_CS3,0<=pdis_CS4<=6.6*X_CS4,
    s_CS1(:,1)==S0_CS1'+0.95*pch_CS1(:,1)-pdis_CS1(:,1)/0.95,
    s_CS2(:,1)==S0_CS2'+0.95*pch_CS2(:,1)-pdis_CS2(:,1)/0.95,
    s_CS3(:,1)==S0_CS3'+0.95*pch_CS3(:,1)-pdis_CS3(:,1)/0.95,
    s_CS4(:,1)==S0_CS4'+0.95*pch_CS4(:,1)-pdis_CS4(:,1)/0.95,
    s_CS1(:,2:24)==s_CS1(:,1:23)+0.95*pch_CS1(:,2:24)-pdis_CS1(:,2:24)/0.95,
    s_CS2(:,2:24)==s_CS2(:,1:23)+0.95*pch_CS2(:,2:24)-pdis_CS2(:,2:24)/0.95,
    s_CS3(:,2:24)==s_CS3(:,1:23)+0.95*pch_CS3(:,2:24)-pdis_CS3(:,2:24)/0.95,
    s_CS4(:,2:24)==s_CS4(:,1:23)+0.95*pch_CS4(:,2:24)-pdis_CS4(:,2:24)/0.95,
    32*0.15<=s_CS1<=32*0.9,32*0.15<=s_CS2<=32*0.9,32*0.15<=s_CS3<=32*0.9,32*0.15<=s_CS4<=32*0.9
    sum(0.95*pch_CS1-pdis_CS1/0.95,2)==32*0.9-S0_CS1',
    sum(0.95*pch_CS2-pdis_CS2/0.95,2)==32*0.9-S0_CS2',
    sum(0.95*pch_CS3-pdis_CS3/0.95,2)==32*0.9-S0_CS3',
    sum(0.95*pch_CS4-pdis_CS4/0.95,2)==32*0.9-S0_CS4'
    ];%Լ������
% %�г���������
Loadcurve=[0.955391944564747,0.978345604157644,1,0.995019488956258,0.972932005197055,0.970333477695972,0.930489389346037,0.890428757037679,0.902771762667822,0.941966219142486,0.911000433087917,0.862061498484192,0.840190558683413,0.831095712429623,0.756604590731919,0.671719359029883,0.611520138588133,0.582936336076224,0.572542226071893,0.574707665656128,0.587267215244695,0.644218276310091,0.755521870939801,0.884798614118666];
PL_base=[5.704;5.705;5.631;6.518;4.890;5.705;5.847]*1000;%���ɷֲ�
PL=PL_base*Loadcurve;%��������
Pf=sdpvar(7,24);%���߹���
Pf(1,:)=PL(1,:)+Pch(1,:)-Pdis(1,:);Pf(2,:)=PL(2,:);Pf(3,:)=PL(3,:);Pf(4,:)=PL(4,:)+Pch(2,:)-Pdis(2,:);Pf(5,:)=PL(5,:)+Pch(3,:)-Pdis(3,:);Pf(6,:)=PL(6,:);Pf(7,:)=PL(7,:)+Pch(4,:)-Pdis(4,:);%���߹������
Pf_limit=1000*[40,40,40,40,40,40,40]';%���߹�������
Pg=sdpvar(10,24);%�����̷ֶε���
Pg_step=1000*[20,5,3,2,2,2,2,2,2,100]';%��������
Price_DSO=[3:12]'*0.1;%�ֶε��
Lagrant_G_left=sdpvar(10,24);%�����̵����½�
Lagrant_G_right=sdpvar(10,24);%�����̵����Ͻ�
b_Lagrant_G_left=binvar(10,24);%�����̵����½粼������
b_Lagrant_G_right=binvar(10,24);%�����̵����Ͻ粼������
Lagrant_L_left=sdpvar(7,24);%��·�����½�
Lagrant_L_right=sdpvar(7,24);%��·�����Ͻ�
b_Lagrant_L_left=binvar(7,24);%��·�����Ͻ粼������
b_Lagrant_L_right=binvar(7,24);%��·�����½粼������
Lagrant_G=sdpvar(1,24);%ƽ��ڵ���
DLMP=sdpvar(7,24);%DLMP
Ckkt=[sum(Pg)==sum(Pf),
    0<=Pg<=Pg_step*ones(1,24),
    Price_DSO*ones(1,24)-Lagrant_G_left+Lagrant_G_right-ones(10,1)*Lagrant_G==0,
    DLMP==ones(7,1)*Lagrant_G+Lagrant_L_right-Lagrant_L_left,
    Pg<=1E6*b_Lagrant_G_left,0<=Lagrant_G_left<=(1-b_Lagrant_G_left),
    Pg_step*ones(1,24)-Pg<=1E6*b_Lagrant_G_right,0<=Lagrant_G_right<=(1-b_Lagrant_G_right),
    -Pf_limit*ones(1,24)<=Pf<=Pf_limit*ones(1,24),
    Pf+Pf_limit*ones(1,24)<=1E6*b_Lagrant_L_left,0<=Lagrant_L_left<=(1-b_Lagrant_L_left),
    -Pf+Pf_limit*ones(1,24)<=1E6*b_Lagrant_L_right,0<=Lagrant_L_right<=(1-b_Lagrant_L_right)];%�г����������KKT����
Obj=sum(sum((Price_DSO*ones(1,24)).*Pg))-sum(sum(DLMP.*PL))+sum(sum((Pg_step*ones(1,24)).*Lagrant_G_right))+sum(sum((Lagrant_L_left+Lagrant_L_right).*(Pf_limit*ones(1,24))));
%�������
C=[Ccs,Ckkt];%Լ������
ops=sdpsettings('solver','gurobi');
result=solvesdp(C,Obj,ops);
DLMP=double(DLMP);
Lagrant_G=double(Lagrant_G);
Pf=double(Pf);
Pch=double(Pch);
Pdis=double(Pdis);
Pg=double(Pg);
%ʵʱ����
load data_RT_potential
load result_order
Pch1=RT_CS1(96).Pch;Pdis1=RT_CS1(96).Pdis;Smin1=RT_CS1(96).Smin;Smax1=RT_CS1(96).Smax;dS1=RT_CS1(96).dS;%ȡ���һ��ʱ�̵Ŀɵ���Ǳ��
Pch2=RT_CS2(96).Pch;Pdis2=RT_CS2(96).Pdis;Smin2=RT_CS2(96).Smin;Smax2=RT_CS2(96).Smax;dS2=RT_CS2(96).dS;
Pch3=RT_CS3(96).Pch;Pdis3=RT_CS3(96).Pdis;Smin3=RT_CS3(96).Smin;Smax3=RT_CS3(96).Smax;dS3=RT_CS3(96).dS;
Pch4=RT_CS4(96).Pch;Pdis4=RT_CS4(96).Pdis;Smin4=RT_CS4(96).Smin;Smax4=RT_CS4(96).Smax;dS4=RT_CS4(96).dS;
Pchmax=[Pch1;Pch2;Pch3;Pch4];
Pdismax=[Pdis1;Pdis2;Pdis3;Pdis4];
Smin=[Smin1;Smin2;Smin3;Smin4];
Smax=[Smax1;Smax2;Smax3;Smax4];
dS=[dS1;dS2;dS3;dS4];
Link=zeros(24,96);%ʱ�λ������(��ǰ1h����Ϊʵʱ15min)
for i=1:24
    Link(i,4*i-3:4*i)=1;
end
DA_Ps=Pch-Pdis;%��ǰ��������
DA_Ps=DA_Ps*Link;%ʱ�λ���
DA_Pf=Pf*Link;%��ǰ�ܹ���
%���վͶ������
Pch=sdpvar(4,96);%���վ��繦��
Pdis=sdpvar(4,96);%���վ�ŵ繦��
S=sdpvar(4,96);%���վ���崢��SOC
Ccs=[0<=Pch<=Pchmax,0<=Pdis<=Pdismax,Smin<=S<=Smax,S(:,1)==0.95*0.25*Pch(:,1)-0.25*Pdis(:,1)/0.95+dS(:,1),S(:,2:96)==S(:,1:95)+0.25*0.95*Pch(:,2:96)-0.25*Pdis(:,2:96)/0.95+dS(:,2:96)];%Լ������
% %�г���������
Loadcurve=[0.955391944564747,0.978345604157644,1,0.995019488956258,0.972932005197055,0.970333477695972,0.930489389346037,0.890428757037679,0.902771762667822,0.941966219142486,0.911000433087917,0.862061498484192,0.840190558683413,0.831095712429623,0.756604590731919,0.671719359029883,0.611520138588133,0.582936336076224,0.572542226071893,0.574707665656128,0.587267215244695,0.644218276310091,0.755521870939801,0.884798614118666];
Loadcurve=Loadcurve*Link;%����96��ʱ��
PL_base=[5.704;5.705;5.631;6.518;4.890;5.705;5.847]*1000;%���ɷֲ�
PL=PL_base*Loadcurve;%��������
Pf=sdpvar(7,96);%���߹���
Pf(1,:)=PL(1,:)+Pch(1,:)-Pdis(1,:);Pf(2,:)=PL(2,:);Pf(3,:)=PL(3,:);Pf(4,:)=PL(4,:)+Pch(2,:)-Pdis(2,:);Pf(5,:)=PL(5,:)+Pch(3,:)-Pdis(3,:);Pf(6,:)=PL(6,:);Pf(7,:)=PL(7,:)+Pch(4,:)-Pdis(4,:);%���߹������
Pf_limit=1000*[40,40,40,40,40,40,40]';%���߹�������
Pg=sdpvar(10,96);%�����̷ֶε���
Pg_step=1000*[20,5,3,2,2,2,2,2,2,100]';%��������
Price_DSO=[3:12]'*0.1;%�ֶε��
Lagrant_G_left=sdpvar(10,96);%�����̵����½�
Lagrant_G_right=sdpvar(10,96);%�����̵����Ͻ�
b_Lagrant_G_left=binvar(10,96);%�����̵����½粼������
b_Lagrant_G_right=binvar(10,96);%�����̵����Ͻ粼������
Lagrant_L_left=sdpvar(7,96);%��·�����½�
Lagrant_L_right=sdpvar(7,96);%��·�����Ͻ�
b_Lagrant_L_left=binvar(7,96);%��·�����Ͻ粼������
b_Lagrant_L_right=binvar(7,96);%��·�����½粼������
Lagrant_G=sdpvar(1,96);%ƽ��ڵ���
DLMP=sdpvar(7,96);%DLMP
Ckkt=[sum(Pg)==sum(Pf),
    0<=Pg<=Pg_step*ones(1,96),
    0.25*Price_DSO*ones(1,96)-Lagrant_G_left+Lagrant_G_right-ones(10,1)*Lagrant_G==0,
    DLMP==ones(7,1)*Lagrant_G+Lagrant_L_right-Lagrant_L_left,
    Pg<=1E6*b_Lagrant_G_left,0<=Lagrant_G_left<=(1-b_Lagrant_G_left),
    Pg_step*ones(1,96)-Pg<=1E6*b_Lagrant_G_right,0<=Lagrant_G_right<=(1-b_Lagrant_G_right),
    -Pf_limit*ones(1,96)<=Pf<=Pf_limit*ones(1,96),
    Pf+Pf_limit*ones(1,96)<=1E6*b_Lagrant_L_left,0<=Lagrant_L_left<=(1-b_Lagrant_L_left),
    -Pf+Pf_limit*ones(1,96)<=1E6*b_Lagrant_L_right,0<=Lagrant_L_right<=(1-b_Lagrant_L_right)];%�г����������KKT����
Obj=0.25*sum(sum((Price_DSO*ones(1,96)).*Pg))-sum(sum(DLMP.*PL))+sum(sum((Pg_step*ones(1,96)).*Lagrant_G_right))+sum(sum((Lagrant_L_left+Lagrant_L_right).*(Pf_limit*ones(1,96))))+0.1*0.25*sum(sum(abs(Pch-Pdis-DA_Ps)));
%�������
C=[Ccs,Ckkt];%Լ������
result=solvesdp(C,Obj);
DLMP=double(DLMP);
Lagrant_G=double(Lagrant_G);
Pf=double(Pf);
Pch=double(Pch);
Pdis=double(Pdis);
Pg=double(Pg);
S=double(S);
double(Obj)
figure(1)%ʵʱ��ǰ����Ա�
hold on
plot(sum(Pf),'r.-')
plot(sum(DA_Pf),'b--');
Cost=sum(sum((Pf-PL).*DLMP));%���õ����
Cost_b=sum(sum(((Pf-PL)>0).*(Pf-PL).*DLMP))%�ܹ������
Cost_s=-sum(sum(((Pf-PL)<0).*(Pf-PL).*DLMP))%���۵�����
Cost_r=0.1*0.25*sum(sum(abs(Pch-Pdis-DA_Ps)))%ƽ�����
MAPE=100*sum(abs((sum(DA_Pf)-sum(Pf))./sum(DA_Pf))/96);%ƽ���ٷ����
%% ���߽ڣ���ǰ�ɵ���Ǳ��Ԥ���㷨�Ա�
%�㷨2��Ⱥ��
clear
clc
close all
load data_EV
X_CS1=data_CS1.X;X_CS2=data_CS2.X;X_CS3=data_CS3.X;X_CS4=data_CS4.X;
S0_CS1=data_CS1.S0;S0_CS2=data_CS2.S0;S0_CS3=data_CS3.S0;S0_CS4=data_CS4.S0;
m1=length(S0_CS1);m2=length(S0_CS2);m3=length(S0_CS3);m4=length(S0_CS4);
X=[sum(X_CS1);sum(X_CS2);sum(X_CS3);sum(X_CS4)];
Energy=[32*0.9*m1-sum(S0_CS1);32*0.9*m2-sum(S0_CS2);32*0.9*m3-sum(S0_CS3);32*0.9*m4-sum(S0_CS4)];
%���վͶ������
Pch=sdpvar(4,24);%���վ��繦��
Pdis=sdpvar(4,24);%���վ�ŵ繦��
Ccs=[0<=Pch<=6.6*X,
    0<=Pdis<=6.6*X,
    sum(Pch-Pdis,2)==Energy
    ];%Լ������
% %�г���������
Loadcurve=[0.955391944564747,0.978345604157644,1,0.995019488956258,0.972932005197055,0.970333477695972,0.930489389346037,0.890428757037679,0.902771762667822,0.941966219142486,0.911000433087917,0.862061498484192,0.840190558683413,0.831095712429623,0.756604590731919,0.671719359029883,0.611520138588133,0.582936336076224,0.572542226071893,0.574707665656128,0.587267215244695,0.644218276310091,0.755521870939801,0.884798614118666];
PL_base=[5.704;5.705;5.631;6.518;4.890;5.705;5.847]*1000;%���ɷֲ�
PL=PL_base*Loadcurve;%��������
Pf=sdpvar(7,24);%���߹���
Pf(1,:)=PL(1,:)+Pch(1,:)-Pdis(1,:);Pf(2,:)=PL(2,:);Pf(3,:)=PL(3,:);Pf(4,:)=PL(4,:)+Pch(2,:)-Pdis(2,:);Pf(5,:)=PL(5,:)+Pch(3,:)-Pdis(3,:);Pf(6,:)=PL(6,:);Pf(7,:)=PL(7,:)+Pch(4,:)-Pdis(4,:);%���߹������
Pf_limit=1000*[40,40,40,40,40,40,40]';%���߹�������
Pg=sdpvar(10,24);%�����̷ֶε���
Pg_step=1000*[20,5,3,2,2,2,2,2,2,100]';%��������
Price_DSO=[3:12]'*0.1;%�ֶε��
Lagrant_G_left=sdpvar(10,24);%�����̵����½�
Lagrant_G_right=sdpvar(10,24);%�����̵����Ͻ�
b_Lagrant_G_left=binvar(10,24);%�����̵����½粼������
b_Lagrant_G_right=binvar(10,24);%�����̵����Ͻ粼������
Lagrant_L_left=sdpvar(7,24);%��·�����½�
Lagrant_L_right=sdpvar(7,24);%��·�����Ͻ�
b_Lagrant_L_left=binvar(7,24);%��·�����Ͻ粼������
b_Lagrant_L_right=binvar(7,24);%��·�����½粼������
Lagrant_G=sdpvar(1,24);%ƽ��ڵ���
DLMP=sdpvar(7,24);%DLMP
Ckkt=[sum(Pg)==sum(Pf),
    0<=Pg<=Pg_step*ones(1,24),
    Price_DSO*ones(1,24)-Lagrant_G_left+Lagrant_G_right-ones(10,1)*Lagrant_G==0,
    DLMP==ones(7,1)*Lagrant_G+Lagrant_L_right-Lagrant_L_left,
    Pg<=1E6*b_Lagrant_G_left,0<=Lagrant_G_left<=(1-b_Lagrant_G_left),
    Pg_step*ones(1,24)-Pg<=1E6*b_Lagrant_G_right,0<=Lagrant_G_right<=(1-b_Lagrant_G_right),
    -Pf_limit*ones(1,24)<=Pf<=Pf_limit*ones(1,24),
    Pf+Pf_limit*ones(1,24)<=1E6*b_Lagrant_L_left,0<=Lagrant_L_left<=(1-b_Lagrant_L_left),
    -Pf+Pf_limit*ones(1,24)<=1E6*b_Lagrant_L_right,0<=Lagrant_L_right<=(1-b_Lagrant_L_right)];%�г����������KKT����
Obj=sum(sum((Price_DSO*ones(1,24)).*Pg))-sum(sum(DLMP.*PL))+sum(sum((Pg_step*ones(1,24)).*Lagrant_G_right))+sum(sum((Lagrant_L_left+Lagrant_L_right).*(Pf_limit*ones(1,24))));
%�������
C=[Ccs,Ckkt];%Լ������
ops=sdpsettings('solver','gurobi');
result=solvesdp(C,Obj,ops);
DLMP=double(DLMP);
Lagrant_G=double(Lagrant_G);
Pf=double(Pf);
Pch=double(Pch);
Pdis=double(Pdis);
Pg=double(Pg);
%ʵʱ����
load data_RT_potential
load result_order
Pch1=RT_CS1(96).Pch;Pdis1=RT_CS1(96).Pdis;Smin1=RT_CS1(96).Smin;Smax1=RT_CS1(96).Smax;dS1=RT_CS1(96).dS;%ȡ���һ��ʱ�̵Ŀɵ���Ǳ��
Pch2=RT_CS2(96).Pch;Pdis2=RT_CS2(96).Pdis;Smin2=RT_CS2(96).Smin;Smax2=RT_CS2(96).Smax;dS2=RT_CS2(96).dS;
Pch3=RT_CS3(96).Pch;Pdis3=RT_CS3(96).Pdis;Smin3=RT_CS3(96).Smin;Smax3=RT_CS3(96).Smax;dS3=RT_CS3(96).dS;
Pch4=RT_CS4(96).Pch;Pdis4=RT_CS4(96).Pdis;Smin4=RT_CS4(96).Smin;Smax4=RT_CS4(96).Smax;dS4=RT_CS4(96).dS;
Pchmax=[Pch1;Pch2;Pch3;Pch4];
Pdismax=[Pdis1;Pdis2;Pdis3;Pdis4];
Smin=[Smin1;Smin2;Smin3;Smin4];
Smax=[Smax1;Smax2;Smax3;Smax4];
dS=[dS1;dS2;dS3;dS4];
Link=zeros(24,96);%ʱ�λ������(��ǰ1h����Ϊʵʱ15min)
for i=1:24
    Link(i,4*i-3:4*i)=1;
end
DA_Ps=Pch-Pdis;%��ǰ��������
DA_Ps=DA_Ps*Link;%ʱ�λ���
DA_Pf=Pf*Link;%��ǰ�ܹ���
%���վͶ������
Pch=sdpvar(4,96);%���վ��繦��
Pdis=sdpvar(4,96);%���վ�ŵ繦��
S=sdpvar(4,96);%���վ���崢��SOC
Ccs=[0<=Pch<=Pchmax,0<=Pdis<=Pdismax,Smin<=S<=Smax,S(:,1)==0.95*0.25*Pch(:,1)-0.25*Pdis(:,1)/0.95+dS(:,1),S(:,2:96)==S(:,1:95)+0.25*0.95*Pch(:,2:96)-0.25*Pdis(:,2:96)/0.95+dS(:,2:96)];%Լ������
% %�г���������
Loadcurve=[0.955391944564747,0.978345604157644,1,0.995019488956258,0.972932005197055,0.970333477695972,0.930489389346037,0.890428757037679,0.902771762667822,0.941966219142486,0.911000433087917,0.862061498484192,0.840190558683413,0.831095712429623,0.756604590731919,0.671719359029883,0.611520138588133,0.582936336076224,0.572542226071893,0.574707665656128,0.587267215244695,0.644218276310091,0.755521870939801,0.884798614118666];
Loadcurve=Loadcurve*Link;%����96��ʱ��
PL_base=[5.704;5.705;5.631;6.518;4.890;5.705;5.847]*1000;%���ɷֲ�
PL=PL_base*Loadcurve;%��������
Pf=sdpvar(7,96);%���߹���
Pf(1,:)=PL(1,:)+Pch(1,:)-Pdis(1,:);Pf(2,:)=PL(2,:);Pf(3,:)=PL(3,:);Pf(4,:)=PL(4,:)+Pch(2,:)-Pdis(2,:);Pf(5,:)=PL(5,:)+Pch(3,:)-Pdis(3,:);Pf(6,:)=PL(6,:);Pf(7,:)=PL(7,:)+Pch(4,:)-Pdis(4,:);%���߹������
Pf_limit=1000*[40,40,40,40,40,40,40]';%���߹�������
Pg=sdpvar(10,96);%�����̷ֶε���
Pg_step=1000*[20,5,3,2,2,2,2,2,2,100]';%��������
Price_DSO=[3:12]'*0.1;%�ֶε��
Lagrant_G_left=sdpvar(10,96);%�����̵����½�
Lagrant_G_right=sdpvar(10,96);%�����̵����Ͻ�
b_Lagrant_G_left=binvar(10,96);%�����̵����½粼������
b_Lagrant_G_right=binvar(10,96);%�����̵����Ͻ粼������
Lagrant_L_left=sdpvar(7,96);%��·�����½�
Lagrant_L_right=sdpvar(7,96);%��·�����Ͻ�
b_Lagrant_L_left=binvar(7,96);%��·�����Ͻ粼������
b_Lagrant_L_right=binvar(7,96);%��·�����½粼������
Lagrant_G=sdpvar(1,96);%ƽ��ڵ���
DLMP=sdpvar(7,96);%DLMP
Ckkt=[sum(Pg)==sum(Pf),
    0<=Pg<=Pg_step*ones(1,96),
    0.25*Price_DSO*ones(1,96)-Lagrant_G_left+Lagrant_G_right-ones(10,1)*Lagrant_G==0,
    DLMP==ones(7,1)*Lagrant_G+Lagrant_L_right-Lagrant_L_left,
    Pg<=1E6*b_Lagrant_G_left,0<=Lagrant_G_left<=(1-b_Lagrant_G_left),
    Pg_step*ones(1,96)-Pg<=1E6*b_Lagrant_G_right,0<=Lagrant_G_right<=(1-b_Lagrant_G_right),
    -Pf_limit*ones(1,96)<=Pf<=Pf_limit*ones(1,96),
    Pf+Pf_limit*ones(1,96)<=1E6*b_Lagrant_L_left,0<=Lagrant_L_left<=(1-b_Lagrant_L_left),
    -Pf+Pf_limit*ones(1,96)<=1E6*b_Lagrant_L_right,0<=Lagrant_L_right<=(1-b_Lagrant_L_right)];%�г����������KKT����
Obj=0.25*sum(sum((Price_DSO*ones(1,96)).*Pg))-sum(sum(DLMP.*PL))+sum(sum((Pg_step*ones(1,96)).*Lagrant_G_right))+sum(sum((Lagrant_L_left+Lagrant_L_right).*(Pf_limit*ones(1,96))))+0.1*0.25*sum(sum(abs(Pch-Pdis-DA_Ps)));
%�������
C=[Ccs,Ckkt];%Լ������
result=solvesdp(C,Obj);
DLMP=double(DLMP);
Lagrant_G=double(Lagrant_G);
Pf=double(Pf);
Pch=double(Pch);
Pdis=double(Pdis);
Pg=double(Pg);
S=double(S);
double(Obj)
figure(1)%ʵʱ��ǰ����Ա�
hold on
plot(sum(Pf),'r.-')
plot(sum(DA_Pf),'b--');
Cost=sum(sum((Pf-PL).*DLMP));%���õ����
Cost_b=sum(sum(((Pf-PL)>0).*(Pf-PL).*DLMP))%�ܹ������
Cost_s=-sum(sum(((Pf-PL)<0).*(Pf-PL).*DLMP))%���۵�����
Cost_r=0.1*0.25*sum(sum(abs(Pch-Pdis-DA_Ps)))%ƽ�����
MAPE=100*sum(abs((sum(DA_Pf)-sum(Pf))./sum(DA_Pf))/96);%ƽ���ٷ����
