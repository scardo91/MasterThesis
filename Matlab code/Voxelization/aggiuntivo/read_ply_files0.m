%read_ply_files(str1,str2)
%
%read PLY files generated from bundler or PMVS into files 




function[coord]=read_ply_files0(str1)

%read the input files
Npts=0;
fp=fopen(str1,'r');
dumstr=fscanf(fp,'%s ',1);
while strcmp(dumstr,'vertex')~=1
    dumstr=fscanf(fp,'%s ',1);
end;
dumstr=fscanf(fp,'%s ',1);
Npts=str2num(dumstr);

dumstr=fscanf(fp,'%s ',1);
nprop=0;
while strcmp(dumstr,'property')
    dumstr=fscanf(fp,'%s ',1);
    dumstr=fscanf(fp,'%s ',1);
    dumstr=fscanf(fp,'%s ',1);
    nprop=nprop+1;
end;



while strcmp(dumstr,'end_header')~=1
    dumstr=fscanf(fp,'%s ',1);
end;


coord=zeros(9,Npts);
for p=1:Npts
    if nprop==10
        ppp=fread(fp,[ 6 1 ],'float');
        coord(1:6,p)=ppp;
        ppp=fread(fp,[ 4 1 ],'uint8');
        coord(7:9,p)=ppp(1:3);
    elseif nprop==9
        ppp=fread(fp,[ 6 1 ],'float');
        coord(1:6,p)=ppp;
        ppp=fread(fp,[ 3 1 ],'uint8');
        coord(7:9,p)=ppp(1:3);
    elseif nprop==7
        ppp=fread(fp,[ 3 1 ],'float');
        coord(1:3,p)=ppp;
        ppp=fread(fp,[ 4 1 ],'uint8');
        coord(7:9,p)=ppp(1:3);
    elseif nprop==6
        ppp=fread(fp,[ 3 1 ],'float');
        coord(1:3,p)=ppp;
        ppp=fread(fp,[ 3 1 ],'uint8');
        coord(7:9,p)=ppp(1:3);
    elseif nprop==3
        ppp=fread(fp,[ 3 1 ],'float');
        coord(1:3,p)=ppp;
        coord(7:9,p)=128;
    end;
end;
fclose(fp);

