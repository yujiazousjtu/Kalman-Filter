%x方向上为匀加速运动，y方向上为匀速运动
%init
initPosCov=100;posMeasErr=20;posModelErr=0.2;
    
x=[0;0;0;0];A=[1,0,0,0;0,1,0,0;0,0,1,0;0,0,0,1];
P=A*5;Q=A*0.4;R=[1,0;0,1]*10;H=[1,0,0,0;0,1,0,0];
K=[0,0;0,0;0,0;0,0];z=[0;0];B=[0,0;0,0;0,0;0,0];
x_=[0;0;0;0];P_=[0,0,0,0;0,0,0,0;0,0,0,0;0,0,0,0];
P=A*initPosCov;R=[1,0;0,1]*posMeasErr;
for k=1:1000
    %setTransitionMat
    dt=1;
    A=[1,0,0,0;0,1,0,0;0,0,1,0;0,0,0,1];
    A(1,3)=dt;A(2,4)=dt;
    sigma_2=3*posModelErr/(dt^3);
    dt3=dt^3*sigma_2/3;dt2=dt^2*sigma_2/2;dt1=dt*sigma_2;  
    for i=1:4
        if i<3
            Q(i,i)=dt3;
        else
            Q(i,i)=dt1;
        end
        if i<3
            Q(i,i+2)=dt2;
            Q(i+2,i)=dt2;
        end
    end
    %predict
    u=[0;0];
    x_=A*x+B*u;P_=A*P*A'+Q;
    %update
    zx(k)=0.5*2*k^2;zy(k)=2*k;
    vx(k)=2*k;vy(k)=2;
    z=[zx(k);zy(k)];
    K=P_*H'*(inv(H*P_*H'+R));
    x=x_+K*(z-H*x_);
    P=P_-K*H*P_;
    %result
    pos=[x(1,1);x(2,1)];pos_x(k)=x(1,1);pos_y(k)=x(2,1);
    vel=[x(3,1);x(4,1)];vel_x(k)=x(3,1);vel_y(k)=x(4,1);
    %disp(vx(k));disp(vel_x(k));
    diff_pos_x(k)=(zx(k)-pos_x(k))/zx(k);
    diff_pos_y(k)=(zy(k)-pos_y(k))/zy(k);
    diff_vel_x(k)=(vx(k)-vel_x(k))/vx(k);
    diff_vel_y(k)=(vy(k)-vel_y(k))/vy(k);
end
%output
figure();
plot(pos_x,'r');
hold on
plot(zx,'b')
legend('x方向滤波位移','x方向理论位移')

figure();
plot(vel_x,'r');
hold on
plot(vx,'b')
legend('x方向滤波速度','x方向理论速度')

figure();
plot(diff_pos_x,'r');
hold on
plot(diff_vel_x,'b')
legend('x方向位移误差','x方向速度误差')

mean1=mean(diff_pos_x);
mean2=mean(diff_vel_x);
disp(mean1);
disp(mean2);

if mean1>=0
    elem1=max(diff_pos_x);
else
    elem1=min(diff_pos_x);
end
if mean2>=0
    elem2=max(diff_vel_x);
else
    elem2=min(diff_vel_x);
end
disp(elem1);
disp(elem2);

figure();
plot(pos_y,'r');
hold on
plot(zy,'b')
legend('y方向滤波位移','y方向理论位移')

figure();
plot(vel_y,'r');
hold on
plot(vy,'b')
legend('y方向滤波速度','y方向理论速度')

figure();
plot(diff_pos_y,'r');
hold on
plot(diff_vel_y,'b')
legend('y方向位移误差','y方向速度误差')

mean3=mean(diff_pos_y);
mean4=mean(diff_vel_y);
disp(mean3);
disp(mean4);

if mean3>=0
    elem3=max(diff_pos_y);
else
    elem3=min(diff_pos_y);
end
if mean4>=0
    elem4=max(diff_vel_y);
else
    elem4=min(diff_vel_y);
end
disp(elem3);
disp(elem4);
