function [x_dec, dec_metric] = polar_dec_BSC(r,K,Rseq),

N = length(r);

n = floor(log2(N));

if (K > N || n ~= log2(N)),
    disp('wrong inputs to function polar_encode');
    return;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Load the 5G standard reliability sequence for 1024
% rel_seq1024;
% 
% %Extract reliability sequence for N
% [ir,jr] = find(Rseq1024 <= N);
% 
% Rseq = Rseq1024(jr);
% 
% %Build polar transform matrix for n
% G = build_G(n);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

max_depth = n+1;

total_nodes = 2^(n+1)-1;

downward_beliefs = zeros(max_depth,N)-100;

downward_beliefs(1,:) = r;

upward_beliefs = zeros(max_depth,N)-100;

message_mat = zeros(total_nodes,total_nodes);

current_node = 1;

while current_node ~= 0, %(sum(upward_beliefs(max_depth,:)== -100) > 0),

    %current_node = current_node

    current_depth = floor(log2(current_node))+1;

    %If the current node is a leaf node
    if current_depth == n+1,
        %pass the belief upward to the parent
        parent_node = floor(current_node/2);

        current_col = current_node - (2^n-1);

        %if the bit is frozen:
        if(sum(Rseq(1:N-K) == current_col) > 0),

            upward_beliefs(current_depth,current_col) = 0;

            message_mat(current_node,parent_node) = 1;

            current_node = parent_node;

        else

            %upward_beliefs(current_depth,current_col) = round(downward_beliefs(current_depth,current_col) < 0);

            upward_beliefs(current_depth,current_col) = round(downward_beliefs(current_depth,current_col) > 0.5);

            message_mat(current_node,parent_node) = 1;

            current_node = parent_node;

        end

    end

    %If the current node is not a leaf node:
    if current_depth < n+1,

        left_child_node = 2*current_node;

        right_child_node = 2*current_node + 1;

        parent_node = floor(current_node/2);

        current_pos = current_node - (2^(current_depth-1)-1);

        %current node's incoming beliefs
        cnvec = (current_pos-1)*(2^(n-current_depth+1))+1:current_pos*(2^(n-current_depth+1));

        cnb = downward_beliefs(current_depth,cnvec);

        lcnvec = length(cnvec);

        %If no message is passed to the left child yet:

        if message_mat(current_node,left_child_node) == 0,

            %downward_beliefs(current_depth+1,cnvec(1:lcnvec/2)) = f(cnb);

            cnb1 = cnb(1:lcnvec/2);

            cnb2 = cnb(lcnvec/2+1:lcnvec);

            %f = min(abs(cnb1),abs(cnb2)).*sign(cnb1).*sign(cnb2);
            
            f = cnb1 .* (1-cnb2) + cnb2 .* (1-cnb1);

            downward_beliefs(current_depth+1,cnvec(1:lcnvec/2)) = f;

            message_mat(current_node,left_child_node) = 1;

            current_node = left_child_node;

        elseif message_mat(current_node,right_child_node) == 0,

            %upward beliefs receieved from left child:

            uwb_lc = upward_beliefs(current_depth+1,cnvec(1:lcnvec/2));

            %downward_beliefs(current_depth+1,cnvec(lcnvec/2+1:lcnvec)) = g(cnb,uwb_lc);

            cnb1 = cnb(1:lcnvec/2);

            cnb2 = cnb(lcnvec/2+1:lcnvec);
            
            fup = uwb_lc .* (1-cnb1) + (1-uwb_lc) .* cnb1;

            %g = fup .* (1-cnb2) + (1-fup) .* cnb2;

            g = (fup .*cnb2)./(fup .*cnb2 + (1-fup) .* (1-cnb2)) ;

            downward_beliefs(current_depth+1,cnvec(lcnvec/2+1:lcnvec)) = g;

            message_mat(current_node,right_child_node) = 1;

            current_node = right_child_node;

        else

            uwb_lc = upward_beliefs(current_depth+1,cnvec(1:lcnvec/2));

            uwb_rc = upward_beliefs(current_depth+1,cnvec(lcnvec/2+1:lcnvec));

            uwb_par = [rem(uwb_lc+uwb_rc,2), uwb_rc];

            upward_beliefs(current_depth,cnvec) = uwb_par;

            if parent_node ~= 0,
                message_mat(current_node,parent_node) = 1;
            end

            current_node = parent_node;

        end

    end

end

u_dec = upward_beliefs(n+1,:);

x_dec = u_dec(Rseq(N-K+1:N));

%dec_metric = sum((downward_beliefs(n+1,Rseq(1:N-K))));

%test_vec = downward_beliefs(n+1,Rseq(1:N-K))

dec_metric = -sum(log(downward_beliefs(n+1,Rseq(1:N-K))+1e-6));


return;