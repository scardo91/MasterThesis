function[]=make_ply_coord(coord,str)

Np=size(coord,2);

ind=setdiff(1:Np,find(isnan(coord(1,:))));
ind=setdiff(ind,find(isnan(coord(2,:))));
ind=setdiff(ind,find(isnan(coord(3,:))));
ind=setdiff(ind,find(isnan(coord(4,:))));
ind=setdiff(ind,find(isnan(coord(5,:))));
ind=setdiff(ind,find(isnan(coord(6,:))));

Np=length(ind);

ffw=fopen(str,'w');
    
fprintf(ffw,'ply\n');
fprintf(ffw,'format ascii 1.0\n');
fprintf(ffw,'element vertex %d\n',Np);
fprintf(ffw,'property float x\n');
fprintf(ffw,'property float y\n');
fprintf(ffw,'property float z\n');
fprintf(ffw,'property float nx\n');
fprintf(ffw,'property float ny\n');
fprintf(ffw,'property float nz\n');
fprintf(ffw,'property uchar diffuse_red\n');
fprintf(ffw,'property uchar diffuse_green\n');
fprintf(ffw,'property uchar diffuse_blue\n');
fprintf(ffw,'end_header\n');

for c=ind
    fprintf(ffw,'%f %f %f 0.0 0.0 0.0 %d %d %d\n',...
			coord(1,c),coord(2,c),coord(3,c),coord(4,c),coord(5,c),coord(6,c));
end;

fclose(ffw);
