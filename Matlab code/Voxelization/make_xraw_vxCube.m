function[]=make_xraw_vxCube(str,vol,col,ncol);

dim=size(vol);
[vvv,vvv_a,vvv_c]=unique(col,'rows');
%vvv=col;
Ncol=size(vvv,1);

if (Ncol>ncol)
    %iii_ccc=floor(colormap(jet(ncol))*255);
    iv=randperm(Ncol);
    iii_ccc=vvv(iv(1:ncol),:);
    iii=zeros(size(vvv,1),1);
    for c=1:size(vvv,1)
        [d,id]=min(max(abs(ones(ncol,1)*vvv(c,:)-iii_ccc),[],2));
        iii(c)=id;
    end;
else
    iii_ccc=vvv;
    [~,iii] = ismember(col,iii_ccc,'rows');
end;
Ncol=min(ncol,size(vvv,1));
iii_ccc=round(iii_ccc);

fp=fopen(str,'wb');
fwrite(fp,[ 'XRAW' ],'uint8');
fwrite(fp,[ 0 4 8 16 ],'uint8');
fwrite(fp,dim,'uint32');
fwrite(fp,Ncol,'uint32');
ind_valid=find(vol(:)>0);
prev=1;
iii=iii(vvv_c);
for c=1:length(ind_valid)
    set_null=prev:ind_valid(c)-1;
    if ~isempty(set_null)
        fwrite(fp,65535*ones(1,length(set_null)),'uint16');
    end;
    fwrite(fp,iii(c)-1,'uint16');
    prev=ind_valid(c)+1;
end;
set_null=prev:dim(1)*dim(2)*dim(3);
if ~isempty(set_null)
        fwrite(fp,65535*ones(1,length(set_null)),'uint16');
end;
%palette
for c=1:size(iii_ccc,1)
    fwrite(fp,uint8([iii_ccc(c,:) 255  ]),'uint8');
end;
fclose(fp);