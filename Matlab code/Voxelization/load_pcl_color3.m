function[vol,col]=load_pcl_color3(str,Nvoxel,Rm,lim0,lim1)

points=read_ply_files(str);
ivalid=intersect(find(isnan(points(1,:))==0),find(isnan(points(2,:))==0));
ivalid=intersect(ivalid,find(isnan(points(3,:))==0));
points=points(:,ivalid);


vol=zeros(Nvoxel,Nvoxel,Nvoxel);
cnt_vol=zeros(Nvoxel,Nvoxel,Nvoxel);
Npts=size(points,2);

qx=linspace(lim0(1)-0.0001,lim1(1)+0.0001,Nvoxel);
qy=linspace(lim0(2)-0.0001,lim1(2)+0.0001,Nvoxel);
qz=linspace(lim0(3)-0.0001,lim1(3)+0.0001,Nvoxel);

ind_to_pos=zeros(1,Npts);
for p=1:Npts
    ix=max(find(qx< points(1,p)));
    iy=max(find(qy< points(2,p)));
    iz=max(find(qz< points(3,p)));
    ip=(ix-1)+(iy-1)*Nvoxel+(iz-1)*(Nvoxel*Nvoxel)+1;
    if ((ip<=(Nvoxel^3))&(~isempty(ix))&(~isempty(iy))&(~isempty(iz)))
        ind_to_pos(p)=ip;
        vol(ix,iy,iz)=1;
    end;
end;

iii=find(ind_to_pos(:)==0);
ind_to_pos(iii)=[];
Npts=length(ind_to_pos)

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