close all;
clear all;

Nvx=512;
Nvx2=Nvx/2;

vet_ply=dir('Dataset/*.ply');

%create .bin files

vetx=zeros(length(vet_ply));
vety=zeros(length(vet_ply));
vetz=zeros(length(vet_ply));

vetx(11)=90;
vety(11)=90;
vetx(1)=225;
vetx(16)=-90;
vety(16)=30;
vetx(2)=90;
vetx(3)=-90;
vetx(4)=90;
vetx(5)=-90;
vety(6)=10;
vetx(7)=180;
vety(7)=180;
vetz(7)=180;
vetx(8)=90;
vety(8)=90;
vetx(9)=225;
vetx(10)=90;
vety(10)=-90;
vetx(12)=180;
vety(12)=180;
vetx(13)=180;
vety(13)=180;
vetz(13)=180;
vetx(14)=90;
vetx(15)=90;
vety(15)=-90;
vetx(17)=90;
vetx(18)=90;
vetx(19)=180;
vety(19)=180;
vetz(19)=270;
vetx(20)=90;
vetx(21)=-90;
vety(21)=90;
vetx(22)=160;
vety(22)=200;
mx=0; my=0; mz=0;
for np=1:length(vet_ply)
    %vetx(np) = 0;
    %vety(np) = 0;
    %vetz(np) = 0;
    str_ply=sprintf('Dataset/%s',vet_ply(np).name);
    str_bin=strrep(str_ply,'.ply','.bin');
    str_col=strrep(str_ply,'.ply','.col');
    str_xraw=strrep(str_ply,'.ply','.xraw');
    %R=rotx(vetx(np))*roty(vety(np))*rotz(vetz(np));
    %R=vetx(np)*vety(np)*vetz(np);
    [vol,col,mx,my,mz]=load_pcl_color2statcol(str_ply,Nvx,np,mx,my,mz);
    %[vol,col]=load_pcl_color2(str_ply,Nvx,R);
    make_xraw_vxCube(str_xraw,vol,floor(col),4000);
    fp=fopen(str_bin,'w');
    ccc=fwrite(fp,vol,'uint8');
    fclose(fp);
     fp2=fopen(str_col,'w');
    ddd=fwrite(fp2,col,'uint8');
    fclose(fp2);
    disp(np);
    
end;
