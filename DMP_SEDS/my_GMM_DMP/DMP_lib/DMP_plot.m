function options = DMP_plot(para,Demo_gnr)
% DMP的绘图函数
% option为绘图类型，既绘制哪些图

if para.options == 1
    figure(1)
    if para.dimension == 2 
        for i=1:para.Ngr_demo+1
            if i==1
                plot(Demo_gnr{i}(1,:),Demo_gnr{i}(2,:),'r');   % 画图：原始轨迹，红色
                hold on
            else
                plot(Demo_gnr{i}(1,:),Demo_gnr{i}(2,:));   % 画图：生成轨迹，蓝色
                hold on
            end
        end
        xlabel('x/mm');
        ylabel('y/mm');
        legend('demons','DMP');
    end
    if para.dimension == 3
        for i=1:para.Ngr_demo+1
            if i==1
                plot3(Demo_gnr{i}(1,:),Demo_gnr{i}(2,:),Demo_gnr{i}(3,:),'r');   % 画图：原始轨迹，红色
                hold on
            else
                plot3(Demo_gnr{i}(1,:),Demo_gnr{i}(2,:),Demo_gnr{i}(3,:),'b');   % 画图：生成轨迹，蓝色
                hold on
            end
        end
        xlabel('x/mm');
        ylabel('y/mm');
        zlabel('z/mm');
        legend('demons','DMP');
    end
    
end

options = para.options;
end