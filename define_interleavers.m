function [intliv_mat,deintliv_mat] = define_interleavers(L,N),

intliv_mat = [];

deintliv_mat = [];

for i=1:L,

    a = randperm(N);

    [b, c] = sort(a);

    intliv_mat = [intliv_mat; a];

    deintliv_mat = [deintliv_mat; c];

end

return;

