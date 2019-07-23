function[vol,col]=load_pcl_color(str,Nvoxel)

points=read_ply_files0(str);

mean_x=mean(points(1,:));
mean_y=mean(points(2,:));
mean_z=mean(points(3,:));

max_x=max(points(1,:));
min_x=min(points(1,:));
delta_x=max(mean_x-min_x,max_x-mean_x);
max_y=max(points(2,:));
min_y=min(points(2,:));
delta_y=max(mean_y-min_y,max_y-mean_y);
max_z=max(points(3,:));
min_z=min(points(3,:));
delta_z=max(mean_z-min_z,max_z-mean_z);
max_x=mean_x+delta_x;
min_x=mean_x-delta_x;
max_y=mean_y+delta_y;
min_y=mean_y-delta_y;
max_z=mean_z+delta_z;
min_z=mean_z-delta_z;

vol=zeros(Nvoxel,Nvoxel,Nvoxel);
cnt_vol=zeros(Nvoxel,Nvoxel,Nvoxel);
Npts=size(points,2);

qx=linspace(min_x-0.0001,max_x+0.0001,Nvoxel);
qy=linspace(min_y-0.0001,max_y+0.0001,Nvoxel);
qz=linspace(min_z-0.0001,max_z+0.0001,Nvoxel);

ind_to_pos=zeros(1,Npts);
for p=1:Npts
    ix=max(find(qx< points(1,p)));
    iy=max(find(qy< points(2,p)));
    iz=max(find(qz< points(3,p)));
    ip=(ix-1)+(iy-1)*Nvoxel+(iz-1)*(Nvoxel*Nvoxel)+1;
    ind_to_pos(p)=ip;
    vol(ix,iy,iz)=1;
end;

ind_to_pos2=zeros(1,Nvoxel*Nvoxel*Nvoxel);
ind=find(vol(:)>0);
ind_to_pos2(ind)=1:length(ind);

col0=zeros(length(ind),4);

for p=1:Npts
    ip=ind_to_pos2(ind_to_pos(p));
    col0(ip,1:3)=col0(ip,1:3)+points(7:9,p)';
    col0(ip,4)=col0(ip,4)+1;
end;

col=col0(:,1:3)./(col0(:,4)*[1 1 1]);