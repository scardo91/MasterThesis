function[vol, col, mean_x, mean_y, mean_z]=load_pcl_color2statcol(str,Nvoxel,ind, mean_x, mean_y, mean_z)

points=read_ply_files(str);
ivalid=intersect(find(isnan(points(1,:))==0),find(isnan(points(2,:))==0));
ivalid=intersect(ivalid,find(isnan(points(3,:))==0));
points=points(:,ivalid);

if ind==1
    mean_x=mean(points(1,:));
    mean_y=mean(points(2,:));
    mean_z=mean(points(3,:));
end

points(1,:)=points(1,:)-mean_x;
points(2,:)=points(2,:)-mean_y;
points(3,:)=points(3,:)-mean_z;


delta_tot=Nvoxel/2;

Max=0+delta_tot;
Min=0-delta_tot;

vol=zeros(Nvoxel,Nvoxel,Nvoxel);

Npts=size(points,2);

q=linspace(Min-0.0001,Max+0.0001,Nvoxel);


ind_to_pos=zeros(1,Npts);
for p=1:Npts
    ix=find(q< points(1,p), 1, 'last' );
    iy=find(q< points(2,p), 1, 'last' );
    iz=find(q< points(3,p), 1, 'last' );
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
