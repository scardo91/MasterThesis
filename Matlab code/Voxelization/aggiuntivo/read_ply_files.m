%read_ply_files(str1,str2)
%
%read PLY files generated from bundler or PMVS into files 




function[coord]=read_ply_files(str1)

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
        ppp=fscanf(fp,'%f ',[ 6 1 ]);
        coord(1:6,p)=ppp;
        ppp=fscanf(fp,'%d ',[ 4 1 ]);
        coord(7:9,p)=ppp(1:3);
    elseif nprop==9
        ppp=fscanf(fp,'%f ',[ 6 1 ]);
        coord(1:6,p)=ppp;
        ppp=fscanf(fp,'%d ',[ 3 1 ]);
        coord(7:9,p)=ppp(1:3);
    elseif nprop==6
        ppp=fscanf(fp,'%f ',[ 3 1 ]);
        coord(1:3,p)=ppp;
        ppp=fscanf(fp,'%d ',[ 3 1 ]);
        coord(7:9,p)=ppp(1:3);
    end;
end;
fclose(fp);

