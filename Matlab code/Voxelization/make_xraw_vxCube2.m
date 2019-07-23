function[]=make_xraw_vxCube2(str,vol,col);

dim=size(vol);

fp=fopen(str,'wb');
fwrite(fp,[ 'XRAW' ],'uint8');
fwrite(fp,[ 0 4 8 0 ],'uint8');
fwrite(fp,dim,'uint32');
%fwrite(fp,0,'uint32');
ind_valid=find(vol(:)>0);
prev=1;
for c=1:length(ind_valid)
    set_null=prev:ind_valid(c)-1;
    if ~isempty(set_null)
        fwrite(fp,zeros(1,length(set_null)*4),'uint8');
    end;
    fwrite(fp,uint8([col(c,:) 255]),'uint8');
    prev=ind_valid(c)+1;
end;
set_null=prev:dim(1)*dim(2)*dim(3);
if ~isempty(set_null)
        fwrite(fp,zeros(1,length(set_null)*4),'uint8');
end;
fclose(fp);