function[vol,col]=load_xraw(str)

fp=fopen(str,'r');
str0=fread(fp,[1 8],'uint8');
dim=fread(fp,[1 3],'uint32');
ncol=fread(fp,1,'uint32');
vol_array=fread(fp,[dim(1) dim(2) dim(3)],'uint16');
col_array=fread(fp,[ ncol 4 ],'uint8');
fclose(fp);

vol=zeros(dim(1),dim(2),dim(3));
ind=find(vol_array(:)~=65535);
vol(ind)=1;
col=zeros(length(ind),3);
col=col_array(vol_array(ind),:);