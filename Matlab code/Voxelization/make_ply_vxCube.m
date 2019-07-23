function[]=make_ply_vxCube(str,vol,col,resol)

if nargin == 2
    ind=find(vol(:)>0);
    col=zeros(length(ind),3);
    col(:,3)=255;
end;

Nvoxel=size(vol,1);

ind=find(vol(:)>0);

Np=length(ind);

ffw=fopen(str,'w');

    
fprintf(ffw,'ply\n');
fprintf(ffw,'format binary_little_endian 1.0\n');
fprintf(ffw,'comment VCGLIB generated\n');
fprintf(ffw,'element vertex %d\n',Np*8);
fprintf(ffw,'property float x\n');
fprintf(ffw,'property float y\n');
fprintf(ffw,'property float z\n');
fprintf(ffw,'property uchar red\n');
fprintf(ffw,'property uchar green\n');
fprintf(ffw,'property uchar blue\n');
fprintf(ffw,'element face %d\n',Np*6);
fprintf(ffw,'property list uchar uint vertex_indices\n');
fprintf(ffw,'element edge 0\n');
fprintf(ffw,'property int vertex1\n');
fprintf(ffw,'property int vertex2\n');
fprintf(ffw,'property uchar red\n');
fprintf(ffw,'property uchar green\n');
fprintf(ffw,'property uchar blue\n');
fprintf(ffw,'end_header\n');


%write vertices
for p=1:length(ind)
    ixp=mod(ind(p)-1,Nvoxel);
    iyp=mod(floor((ind(p)-1)/Nvoxel),Nvoxel);
    izp=floor((ind(p)-1)/(Nvoxel*Nvoxel));
    
    fwrite(ffw,[ ixp iyp izp ]*resol,'float32');
    fwrite(ffw,col(p,:),'uint8');
    fwrite(ffw,[ (ixp+1) iyp izp ]*resol,'float32');
    fwrite(ffw,col(p,:),'uint8');
    fwrite(ffw,[ ixp (iyp+1) izp ]*resol,'float32');
    fwrite(ffw,col(p,:),'uint8');
    fwrite(ffw,[ (ixp+1) (iyp+1) izp ]*resol,'float32');
    fwrite(ffw,col(p,:),'uint8');
    fwrite(ffw,[ ixp iyp (izp+1) ]*resol,'float32');
    fwrite(ffw,col(p,:),'uint8');
    fwrite(ffw,[ (ixp+1) iyp (izp+1) ]*resol,'float32');
    fwrite(ffw,col(p,:),'uint8');
    fwrite(ffw,[ ixp (iyp+1) (izp+1) ]*resol,'float32');
    fwrite(ffw,col(p,:),'uint8');
    fwrite(ffw,[ (ixp+1) (iyp+1) (izp+1) ]*resol,'float32');
    fwrite(ffw,col(p,:),'uint8');
end;
    
for p=1:length(ind)
    ixp=mod(ind(p)-1,Nvoxel)+1;
    iyp=mod(floor((ind(p)-1)/Nvoxel),Nvoxel)+1;
    izp=floor((ind(p)-1)/(Nvoxel*Nvoxel))+1;
    
    ip=p-1;
    fwrite(ffw,4,'uint8');
    fwrite(ffw,[ (ip*8) (ip*8+1) (ip*8+2) (ip*8+3) ],'uint32');
    fwrite(ffw,4,'uint8');
    fwrite(ffw,[ (ip*8) (ip*8+1) (ip*8+4) (ip*8+5) ],'uint32');
    fwrite(ffw,4,'uint8');
    fwrite(ffw,[ (ip*8) (ip*8+2) (ip*8+4) (ip*8+6) ],'uint32');
    fwrite(ffw,4,'uint8');
    fwrite(ffw,[ (ip*8+1) (ip*8+3) (ip*8+5) (ip*8+7) ],'uint32');
    fwrite(ffw,4,'uint8');
    fwrite(ffw,[ (ip*8+4) (ip*8+5) (ip*8+6) (ip*8+7) ],'uint32');
    fwrite(ffw,4,'uint8');
    fwrite(ffw,[ (ip*8+2) (ip*8+3) (ip*8+6) (ip*8+7) ],'uint32');
end;
    
   
fclose(ffw);
