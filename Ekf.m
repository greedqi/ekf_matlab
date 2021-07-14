% Author:lsq
% Data:7.13
% LastEditTime: 2021.7.14
% Description: EKF滤波类
classdef Ekf
    properties
        X;%状态向量，行向量，所有使用X的地方应该转置。
        P;
        H;
        R;
        Q;
        init_flag=false;%初始化完成后为true
        last_t;
    end
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Name:Init
        % Author:lsq
        % Data:7.13
        % Description :EKF初始化，在Update中调用
        % param:z:（行向量）
        %       t:当前绝对时间单位秒  
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj=Init(obj,z,t)
            if obj.init_flag
                return;
            end
            obj.init_flag=true;
            obj.last_t=t;
            obj.X=z;
            
            obj.H=eye(size(obj.X,2));
            
            obj.P=eye(size(obj.X,2))*1e-08;
            
            %测量噪声协方差矩阵
            obj.R=[ 1.0e-04  0  0  0  0;
                0  1.0e-04  0  0  0;
                0  0  1.0e-06  0  0;
                0  0  0  2.0e-04  0;
                0  0  0  0  1.0e-05;];
            %模型噪声协方差矩阵
            obj.Q=[ 1.0e-3 0   0     0   0;
                0   1.0e-3 0     0   0;
                0   0   1.0e-3 0   0;
                0   0   0     5e-3 0;
                0   0   0     0   2.5e-4;];
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Name:Predict
        % Author:lsq
        % Data:7.13
        % Description :EKF预测，在Update中调用，
        %              此段代码应遵循开闭原则，当控制对象模型改变时，确保只需修改此函数。
        % param:u:控制量（行向量）
        %       t:当前绝对时间单位秒  
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj=Predict(obj,u,t)
            L=4;
            
            dt=t-obj.last_t;
            x=obj.X(1);
            y=obj.X(2);
            phi=obj.X(3);
            v=obj.X(4);
            delta=obj.X(5);
            a=u(1);
            
            p_x=x+cos(phi)*cos(delta)*v*dt;
            p_y=y+sin(phi)*cos(delta)*v*dt;
            p_phi=phi+2*v*sin(delta)*dt/L;
            p_v=v;
            p_delta=delta;
            
            obj.X=[p_x p_y p_phi p_v p_delta];
            
            F =[
                1 0  -(v + a * dt / 2) * sin(phi) * cos(delta) * dt  cos(delta) * cos(phi) * dt -(v + a * dt / 2) * cos(phi) * sin(delta) * dt;
                0 1 (v + a * dt / 2) * cos(phi) * cos(delta) * dt    cos(delta) * sin(phi) * dt -(v + a * dt / 2) * sin(phi) * sin(delta) * dt;
                0 0 1 2 * sin(delta) / L * dt  2 * v * cos(delta) * dt / L;
                0 0 0 1 0;
                0 0 0 0 1; ];
            
            obj.P=F*obj.P*F'+obj.Q;
            
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Name:Update
        % Author:lsq
        % Data:7.13
        % Description :EKF Update
        %              给外部调用ekf的接口。
        % param:z:测量值（行向量）
        %       u:控制量（行向量）
        %       t:当前绝对时间单位秒  
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj=Update(obj,z,u,t)
            
            obj=obj.Init(z,t);
            obj=obj.Predict(u,t);
            K = obj.P * obj.H'/(obj.H * obj.P * obj.H' + obj.R);
            
            y= z' - obj.H * obj.X';%X，z是行向量，转置。
            obj.X= obj.X + (K * y)';
            obj.P=(eye(size(obj.X,2)) - K * obj.H) * obj.P;
        end
        
        function X=Get_x(obj)
            X=obj.X;
        end
        function P=Get_P(obj)
            P=obj.P;
        end
        
    end
end