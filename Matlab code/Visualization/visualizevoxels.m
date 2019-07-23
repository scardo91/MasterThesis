clear all;
%close all;

fp = fopen('tf_color.col','rb');  %
col = fread(fp,'double');
%col = fread(fp,'uint8');
fclose(fp);
col = reshape(col,[3, numel(col)/3]);
col = rot90(col,3);
col = fliplr(col);
%col = double(col);

fp2 = fopen('frame0001.bin','rb');
vol = fread(fp2,'uint8');
fclose(fp2);

sz = 2^(log2(numel(vol))/3);
vol = reshape(vol,[sz sz sz]);
figure();
[X,Y,Z] = ind2sub(size(vol),find(vol));
scatter3(X,Y,Z,1,col(1:numel(X),:)./255);
axis([1 size(vol, 1) 1 size(vol,2) 1 size(vol,3)]);
axis square

%%uncomment if orig color is available
fp = fopen('orig_color.col','rb');
o_col = fread(fp,'double');
fclose(fp);
o_col = reshape(o_col,[3, numel(o_col)/3]);
o_col = rot90(o_col,3);
o_col = fliplr(o_col);
figure()
scatter3(X,Y,Z,1,o_col(1:numel(X),:)./255);
axis([1 size(vol, 1) 1 size(vol,2) 1 size(vol,3)]);

axis square