function[vol]=load_pcl(str,Nvoxel)

points=read_ply_files0(str);
ind=find(~isnan(points(1,:)));
ind=intersect(ind,find(~isnan(points(2,:))));
ind=intersect(ind,find(~isnan(points(3,:))));
points=points(:,ind);

mean_x=mean(points(1,:));
mean_y=mean(points(2,:));
mean_z=mean(points(3,:));

delta_voxel=0;
max_x=max(points(1,:));
min_x=min(points(1,:));
delta_x=max(mean_x-min_x,max_x-mean_x);
if (delta_voxel<delta_x)
    delta_voxel=delta_x;
end;
max_y=max(points(2,:));
min_y=min(points(2,:));
delta_y=max(mean_y-min_y,max_y-mean_y);
if (delta_voxel<delta_y)
    delta_voxel=delta_y;
end;
max_z=max(points(3,:));
min_z=min(points(3,:));
delta_z=max(mean_z-min_z,max_z-mean_z);
if (delta_voxel<delta_z)
    delta_voxel=delta_z;
end;
max_x=mean_x+delta_voxel;
min_x=mean_x-delta_voxel;
max_y=mean_y+delta_voxel;
min_y=mean_y-delta_voxel;
max_z=mean_z+delta_voxel;
min_z=mean_z-delta_voxel;

vol=zeros(Nvoxel,Nvoxel,Nvoxel);
cnt_vol=zeros(Nvoxel,Nvoxel,Nvoxel);
Npts=size(points,2);

qx=linspace(min_x-0.0001,max_x+0.0001,Nvoxel);
qy=linspace(min_y-0.0001,max_y+0.0001,Nvoxel);
qz=linspace(min_z-0.0001,max_z+0.0001,Nvoxel);

for p=1:Npts
    ix=max(find(qx< points(1,p)));
    iy=max(find(qy< points(2,p)));
    iz=max(find(qz< points(3,p)));
    vol(ix,iy,iz)=1;
end;

% 
% figure(1);
% scatter3(points(1,:),points(2,:),points(3,:));