function[bitstream,prob_mat]=octree_coding(vol)

dim_vet=size(vol);

if numel(dim_vet) ~= 3
    error('Input must be a 3D matrix in single or double precision');
end

if ((dim_vet(1)~=dim_vet(2))|(dim_vet(3)~=dim_vet(2)))
    error('Dimensions of the 3D matrix must be the same');
end;

if (mod(log2(dim_vet(1)),1)~=0)
    error('Dimensions of the 3D matrix must be power of two');
end;

Ndecomp=log2(dim_vet(1));

bitstream = code_octet(vol);

end


function bool_empty = is_cube_notempty(c)
    if length(find(c(:)>0))>0
        bool_empty=1;
    else
        bool_empty=0;
    end;
end

function vet_bin = code_octet(v)
% MYMEAN Example of a local function.
dim_oct=size(v,1);
dim_oct_half=dim_oct/2;


cube0=v(1:dim_oct_half,1:dim_oct_half,1:dim_oct_half);
cube1=v(dim_oct_half+1:dim_oct,1:dim_oct_half,1:dim_oct_half);
cube2=v(1:dim_oct_half,dim_oct_half+1:dim_oct,1:dim_oct_half);
cube3=v(dim_oct_half+1:dim_oct,dim_oct_half+1:dim_oct,1:dim_oct_half);
cube4=v(1:dim_oct_half,1:dim_oct_half,dim_oct_half+1:dim_oct);
cube5=v(dim_oct_half+1:dim_oct,1:dim_oct_half,dim_oct_half+1:dim_oct);
cube6=v(1:dim_oct_half,dim_oct_half+1:dim_oct,dim_oct_half+1:dim_oct);
cube7=v(dim_oct_half+1:dim_oct,dim_oct_half+1:dim_oct,dim_oct_half+1:dim_oct);

vet_bin=[ is_cube_notempty(cube0) is_cube_notempty(cube1) is_cube_notempty(cube2)];
vet_bin=[ vet_bin is_cube_notempty(cube3) is_cube_notempty(cube4) is_cube_notempty(cube5)];
vet_bin=[ vet_bin is_cube_notempty(cube6) is_cube_notempty(cube7) ];

if (dim_oct>2)
    vet_bin_add=[];
    for c=0:7
        if vet_bin(c+1)==1
            eval(sprintf('vv=code_octet(cube%d);',c));
            vet_bin_add=[vet_bin_add vv];
        end;
    end;
    vet_bin=[ vet_bin vet_bin_add ];
end;


end