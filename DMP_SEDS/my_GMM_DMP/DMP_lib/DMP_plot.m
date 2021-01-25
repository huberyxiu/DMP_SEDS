function options = DMP_plot(para,Demo_gnr)
% DMP�Ļ�ͼ����
% optionΪ��ͼ���ͣ��Ȼ�����Щͼ

if para.options == 1
    figure(1)
    if para.dimension == 2 
        for i=1:para.Ngr_demo+1
            if i==1
                plot(Demo_gnr{i}(1,:),Demo_gnr{i}(2,:),'r');   % ��ͼ��ԭʼ�켣����ɫ
                hold on
            else
                plot(Demo_gnr{i}(1,:),Demo_gnr{i}(2,:));   % ��ͼ�����ɹ켣����ɫ
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
                plot3(Demo_gnr{i}(1,:),Demo_gnr{i}(2,:),Demo_gnr{i}(3,:),'r');   % ��ͼ��ԭʼ�켣����ɫ
                hold on
            else
                plot3(Demo_gnr{i}(1,:),Demo_gnr{i}(2,:),Demo_gnr{i}(3,:),'b');   % ��ͼ�����ɹ켣����ɫ
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