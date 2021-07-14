close all;clear;clc;
ekf=Ekf();
data=textread('C:\Users\lisiqi\Desktop\ekf.txt');
rows=size(data,1);

sensor_data=data(1:2:rows,:);
ekf_data=data(2:2:rows,:);
ekf_local_data=zeros(rows/2,5);
p_data=zeros(rows/2-1,5);

ekf_rows=size(ekf_data,1);
seonsor_rows=size(sensor_data,1);

% plot(sensor_data(:,1),sensor_data(:,2),'-oy',ekf_data(:,2),ekf_data(:,3),'-oc')
% xlim([-0.6-0.5 1.95-0.6-0.5])
% ylim([-0.1 1.85])

l=4;
%phi转换成相对phi（-inf，+inf）
for i = 2:ekf_rows
    last_phi=sensor_data(i-1,3);
    phi=sensor_data(i,3);
    if phi-last_phi>pi
        phi=phi-2*pi;
    elseif phi-last_phi<-pi
        phi=phi+2*pi;
    end
    sensor_data(i,3)=phi;
end


%ekf滤波
A_P=zeros(5,5);
for i = 1:seonsor_rows
    
    t=ekf_data(i,1);
    x=sensor_data(i,1);
    y=sensor_data(i,2);
    phi=sensor_data(i,3);
    v=sensor_data(i,4);
    delta=sensor_data(i,5);
    
    z=[x y phi v delta];
    a=ekf_data(i,8);
    u=[a];
    
    ekf=ekf.Update(z,u,t/1000);
    ekf_local_data(i,:)=ekf.Get_x();
    p_data(i,:)=ekf.Get_X_();
    A_P=A_P+ekf.Get_P();
    

end
average_P=A_P/seonsor_rows

%计算模型预测的点存到p_data中，与ekf无关debug用
% for i = 1:ekf_rows-1
%
%     dt=(ekf_data(i+1,1)-ekf_data(i,1))/1000;
%
%     x=ekf_local_data(i,1);
%     y=ekf_local_data(i,2);
%     phi=ekf_local_data(i,3);
%     v=ekf_local_data(i,4);
%     delta=ekf_local_data(i,5);
%
%     a=ekf_data(i,8);
%
%     p_x=x+cos(phi)*cos(delta)*v*dt;
%     p_y=y+sin(phi)*cos(delta)*v*dt;
%     p_phi=phi+2*v*sin(delta)*dt/l;
%     p_v=v;
%     p_delta=delta;
%
%     p_data(i,:)=[p_x p_y p_phi p_v p_delta];
%
% end


% %计算ekf模型预测与测量值的协方差。与ekf无关debug用
predict_senor_covariance=zeros(5,5);
for i=1:ekf_rows-1
    for j=1:4
        for k=1:4
            predict_senor_covariance(j,k)=predict_senor_covariance(j,k)+(p_data(i,j)-sensor_data(i,j))*(p_data(i,k)-sensor_data(i,k));
            
        end
    end
    if(p_data(i,1)>3)
        i
    end
    
end
re=predict_senor_covariance/(ekf_rows-1)




%使用静态下的测量矩阵计算静态时的R
% R=zeros(4,4);
% temp=[0.05 0.05 5/180*pi 0];
% for i=1:seonsor_rows
%
%     for j=1:4
%         for k=1:4
%
%             R(j,k)=R(j,k)+(sensor_data(i,j)-temp(j))*(sensor_data(i,k)-temp(k));
%
%         end
%     end
%
% end
% R=R/(seonsor_rows);

hold on
plot(p_data(:,1),p_data(:,2),'xg',sensor_data(:,1),sensor_data(:,2),'-xr',ekf_local_data(:,1),ekf_local_data(:,2),'-xk');



