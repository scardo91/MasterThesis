function[vol]=octree_decoding(bit_str,Ndim)

dim_str=length(bit_str);

[vol , pos]= decode_octet(bit_str,Ndim,0);

end


 function[cube_oct,pos] = decode_octet(str,dim_oct,pos)
% MYMEAN Example of a local function.

dim_oct_half=dim_oct/2;
cube_oct=zeros(dim_oct,dim_oct,dim_oct);
head_str=str(pos+1:pos+8);

if dim_oct==2
    cube_oct(1,1,1)=head_str(1);
    cube_oct(2,1,1)=head_str(2);
    cube_oct(1,2,1)=head_str(3);
    cube_oct(2,2,1)=head_str(4);
    cube_oct(1,1,2)=head_str(5);
    cube_oct(2,1,2)=head_str(6);
    cube_oct(1,2,2)=head_str(7);
    cube_oct(2,2,2)=head_str(8);
    pos=pos+8;
else
    pos=pos+8;
    if (head_str(1)>0)
        [dec_cube,pos0]=decode_octet(str,dim_oct_half,pos);
        cube_oct(1:dim_oct_half,1:dim_oct_half,1:dim_oct_half)=dec_cube;
        pos=pos0;
    end;
    if (head_str(2)>0)
        [dec_cube,pos0]=decode_octet(str,dim_oct_half,pos);
        cube_oct(dim_oct_half+1:dim_oct,1:dim_oct_half,1:dim_oct_half)=dec_cube;
        pos=pos0;
    end;
    if (head_str(3)>0)
        [dec_cube,pos0]=decode_octet(str,dim_oct_half,pos);
        cube_oct(1:dim_oct_half,dim_oct_half+1:dim_oct,1:dim_oct_half)=dec_cube;
        pos=pos0;
    end;
    if (head_str(4)>0)
        [dec_cube,pos0]=decode_octet(str,dim_oct_half,pos);
        cube_oct(dim_oct_half+1:dim_oct,dim_oct_half+1:dim_oct,1:dim_oct_half)=dec_cube;
        pos=pos0;
    end;
    
    if (head_str(5)>0)
        [dec_cube,pos0]=decode_octet(str,dim_oct_half,pos);
        cube_oct(1:dim_oct_half,1:dim_oct_half,dim_oct_half+1:dim_oct)=dec_cube;
        pos=pos0;
    end;
    if (head_str(6)>0)
        [dec_cube,pos0]=decode_octet(str,dim_oct_half,pos);
        cube_oct(dim_oct_half+1:dim_oct,1:dim_oct_half,dim_oct_half+1:dim_oct)=dec_cube;
        pos=pos0;
    end;
    if (head_str(7)>0)
        [dec_cube,pos0]=decode_octet(str,dim_oct_half,pos);
        cube_oct(1:dim_oct_half,dim_oct_half+1:dim_oct,dim_oct_half+1:dim_oct)=dec_cube;
        pos=pos0;
    end;
    if (head_str(8)>0)
        [dec_cube,pos0]=decode_octet(str,dim_oct_half,pos);
        cube_oct(dim_oct_half+1:dim_oct,dim_oct_half+1:dim_oct,dim_oct_half+1:dim_oct)=dec_cube;
        pos=pos0;
    end;
end;




end