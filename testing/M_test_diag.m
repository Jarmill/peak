load('M_test.mat');

%testing block-diagonalization of moment matrices coming from an
%equivariant symmetry x->-x. This symmetry arises in linear systems and the
%sym_attractor. Hopefully this structure may be used on linear systems
%analysis without too much pain. The dual (function) problem would have
%v(x) as an even polynomial.


ind_even = [];
ind_odd = [];
step = 0;
curr = 0;
curr_seq = [];
for d = 0:order
   %number of monomials of degree d
   curr = nchoosek(d+n-1, d);
   curr_seq = [curr_seq curr];
   if  mod(d, 2) == 0
       %even       
       ind_even = [ind_even, step + (1:curr)];
   else
       %odd
       ind_odd = [ind_odd, step + (1:curr)];
   end
   step = step + curr;
end

%occupation measure is of higher degree (by 1?)
d = order + 1;
curr = nchoosek(d+n-1, d);
curr_seq = [curr_seq curr];
if  mod(d, 2) == 0
   %even       
   ind_odd_occ = ind_odd;
   ind_even_occ = [ind_even, step + (1:curr)];
else
   %odd
   ind_odd_occ = [ind_odd, step + (1:curr)];
   ind_even_occ = ind_even;
end
    

M0_even = M0(ind_even, ind_even);
M0_odd = M0(ind_odd, ind_odd);
M0_blk = blkdiag(M0_even, M0_odd);

Mp_even = Mp(ind_even, ind_even);
Mp_odd = Mp(ind_odd, ind_odd);
Mp_blk = blkdiag(Mp_even, Mp_odd);

Mocc_even = Mocc(ind_even_occ, ind_even_occ);
Mocc_odd = Mocc(ind_odd_occ, ind_odd_occ);
Mocc_blk = blkdiag(Mocc_even, Mocc_odd);
%odd is replicated in even