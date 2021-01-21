%修改B矩阵，加入加速度输入
%init
initPosCov=100;posMeasErr=20;posModelErr=0.2;
    
x=[0;0;0;0];A=[1,0,0,0;0,1,0,0;0,0,1,0;0,0,0,1];
P=A*5;Q=A*0.4;R=[1,0;0,1]*10;H=[1,0,0,0;0,1,0,0];
K=[0,0;0,0;0,0;0,0];z=[0;0];B=[0,0;0,0;0,0;0,0];
x_=[0;0;0;0];P_=[0,0,0,0;0,0,0,0;0,0,0,0;0,0,0,0];
P=A*initPosCov;R=[1,0;0,1]*posMeasErr;
for k=1:1000
    %setTransitionMat
    dt=1;ax=10;ay=0;
    A=[1,0,0,0;0,1,0,0;0,0,1,0;0,0,0,1];
    A(1,3)=dt;A(2,4)=dt;
    B=[(k-1)*dt+0.5*dt^2,0;0,(k-1)*dt+0.5*dt^2;dt,0;0,dt];
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
    u=[ax;ay];
    x_=A*x+B*u;P_=A*P*A'+Q;
    %disp(x(3,1));disp(x_(3,1));
    %update
    zx(k)=0.5*ax*k^2;zy(k)=2*k;
    vx(k)=ax*k;vy(k)=2;
    z=[zx(k);zy(k)];
    K=P_*H'*(inv(H*P_*H'+R));
    x=x_+K*(z-H*x_);
    P=P_-K*H*P_;
    %result
    pos=[x(1,1);x(2,1)];pos_x(k)=x(1,1);pos_y(k)=x(2,1);
    %x(3,1)=vx(k);x(4,1)=vy(k);
    vel=[x(3,1);x(4,1)];vel_x(k)=x(3,1);vel_y(k)=x(4,1);
    diff_pos(k)=(zx(k)-pos_x(k))/zx(k);
    diff_vel(k)=(vx(k)-vel_x(k))/vx(k);
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
plot(diff_pos,'r');
hold on
plot(diff_vel,'b')
legend('x方向位移误差','x方向速度误差')

mean1=mean(diff_pos);
mean2=mean(diff_vel);
disp(mean1);
disp(mean2);

min1=min(diff_pos);
max2=max(diff_vel);
disp(min1);
disp(max2);
