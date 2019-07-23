close all;
clear all;

cur=zeros(8,8,8);
cur(1:4,1:4,1:4)=1;
cur(1:2,1:2,5)=1;
ref=zeros(8,8,8);
ref(3:6,3:6,1:4)=1;

voxelPlot(cur);
title('cur');
voxelPlot(ref);
title('ref');

[pred_array,mv]=find_voxel_pred(uint8(cur(:)),uint8(ref(:)),8,4,4);

pred=double(reshape(pred_array,[8 8 8]));

voxelPlot(pred);
title('pred');
hold on;
cnt=1;
for z=1:4:8
    for y=1:4:8
        for x=1:4:8
            mv0=mv(:,cnt);
            x1=x+mv0(1);
            y1=y+mv0(2);
            z1=z+mv0(3);
            line([x x1],[y y1],[z z1],'Color','red','LineWidth',4);
            a=plot3(x1,y1,z1,'b.');
            set(a,'Color','r');
            set(a,'MarkerSize',30);
            cnt=cnt+1;
        end;
    end;
end;
hold off;
