clc;
clear;

%Parameters of the outer code: RS code
m = 6;           % Number of bits per symbol
n = 2^m - 1;     % Codeword length
%k = 225; %225;           % Message length

%Number of information bits per whole block
%lbits = m*k;

%number of segments
seg_num = 6;

%number of bits required to encode the index of each segment
ind_bits_num = ceil(log2(seg_num));

%number of information bits in each slice
k_seg = ceil(n*m/seg_num);

%Outer code : polar codeword length
n_seg = 128;

%Load the 5G standard reliability sequence for 1024
rel_seq1024;

%Extract reliability sequence for n_seg
[ir,jr] = find(Rseq1024 <= n_seg);

Rseq = Rseq1024(jr);

%Build polar transform matrix
n_polar = floor(log2(n_seg));

G = build_G(n_polar);

%This flag bit denotes whether explicit indexing is used
%ind = 0

k_vec = 53:10:53;

p_bsc_vec = 0.025:0.005:0.050;

BER0 = ones(length(k_vec),length(p_bsc_vec));

BER1 = ones(length(k_vec),length(p_bsc_vec));

FER0 = ones(length(k_vec),length(p_bsc_vec));

FER1 = ones(length(k_vec),length(p_bsc_vec));

SER0 = ones(length(k_vec),length(p_bsc_vec));

SER1 = ones(length(k_vec),length(p_bsc_vec));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




for ind = 0:1,

    %AWGN variance
    %SNR_db = 1;

    for ik = 1:length(k_vec),

        k = k_vec(ik);

        lbits = m*k;

        for ipb = 1:length(p_bsc_vec),

            p_bsc = p_bsc_vec(ipb);

            %sigma_noise = sqrt(0.5*10^(-SNR_db/10));

            bit_errors = 0;

            frame_errors = 0;

            symbol_errors = 0;

            sent_blocks = 5000;

            %coset_mat = rand(seg_num , n_seg) < 0.1;
            
            %Hcs = round(hadamard(n_seg) > 0);

            %coset_mat = Hcs(2:seg_num+1,:);

            %ercs = ones(seg_num,1) * round(rand(1,n_seg));

            %coset_mat = rem(Hcs(1:seg_num,:) + ercs, 2);

            %coset_mat = seqs_max_dist(seg_num,n_seg,10,20,0.1);

            %coset_mat = select_coset_leaders(seg_num,80,Rseq);

            
            for sb = 1:sent_blocks,



                %Define interleavers and de-interleavers for matched decoder scheme
                %[intliv_mat,deintliv_mat] = define_interleavers(seg_num,n_seg);

                %Define coset leaders for matched decoder scheme
                
%                 coset_mat = [];
% 
%                 for snco = 0:seg_num-1,
%                     
%                     u_snco = zeros(1,n_seg);
% 
%                     %selected_pos = Rseq(n_seg-k_seg-ind_bits_num+1:n_seg-k_seg);
% 
%                     %u_snco(selected_pos) = gf2bit(snco,ind_bits_num);
% 
%                     selected_pos = Rseq(1:n_seg-k_seg); %-ind_bits_num+1:n_seg-k_seg);
% 
%                     u_snco(selected_pos) = round(rand(size(selected_pos)));
% 
% 
%                     cs_leader = rem(u_snco * G , 2);
% 
%                     coset_mat = [coset_mat ; cs_leader];
% 
%                 end

                coset_mat = rand(seg_num,n_seg) < 0.5;

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    

                %generating information bits
                u = round(rand(1,lbits) < 0.5);

                %Converting bits to m-bit symbols
                unb = bit2gf(u,m);

                %converting symbols to MATLAB's GF(2^m) data format
                msg = gf(unb,m);

                %Encoding by RS code
                code = rsenc(msg,n,k);

                %Converting the codeword from MATLAB's GF(2^m) format to integer format
                codeint = zeros(1,n);

                for j=0:2^m-1,

                    codeint = codeint + j*(code == j);

                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                %Converting the symbols sequences to bit sequence
                x = gf2bit(codeint,m);

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %Slicing the sequence into segments

                %padding the input sequence to get an integer number of segments
                if length(x) < k_seg*seg_num,

                    x = [x, zeros(1,k_seg*seg_num-length(x))];
                end

                %encoding each slice by a polar code

                polar_cws = [];

                for si = 1:seg_num,

                    xi = x(1+(si-1)*k_seg:si*k_seg);

                    if ind == 1,

                        xi = [xi, gf2bit(si-1,ind_bits_num)];
                    end

                    %Calling the polar code function
                    xi_out = polar_encode(xi,n_seg,G,Rseq);
                    %%%%%

                    if ind == 1,
                        polar_cws = [polar_cws; xi_out];
                    end


                    %interleave each segment if matched decoding is employed
                    if ind == 0,

                        %polar_cws = [polar_cws; xi_out(intliv_mat(si,:))];

                        polar_cws = [polar_cws; rem(xi_out+coset_mat(si,:),2)];
                    end


                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %Passing through a BSC channel
                
                z_mat = (rand(size(polar_cws)) < p_bsc);

                y_mat = rem(polar_cws+z_mat, 2);

                %Shuffle the noisy observations
                shuf = randperm(seg_num);

                y_mat = y_mat(shuf,:);

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%% Decoding the segments and detecting the indexes to put them in order
                x_rec_mat = zeros(seg_num,k_seg);

                %decoding in case of explicit indexing
                if ind == 1,
                    for si=1:seg_num,

                        %calling the polar decoder
                        k_val = k_seg+ind_bits_num;
                        
                        %a priori probability of message bits being "1"
                        %f_vec = 0.5*ones(1,n_seg);
                        %a priori reliability of frozen bits being "1" is zero
                        %f_vec(Rseq(1:n_seg-k_seg)) = 0;
                        
                        out_rel_vec = (1-p_bsc)*y_mat(si,:)+(p_bsc)*(1-y_mat(si,:));

                        %xi_rec = polar_dec_BSC(y_mat(si,:),k_seg+ind_bits_num,Rseq);

                        xi_rec = polar_dec_BSC(out_rel_vec,k_val,Rseq);

                        %[out_rec, xi_rec] = polar_decode(out_rel_vec,f_vec);

                        xi_ind_bits = xi_rec(end-ind_bits_num+1:end);

                        xi_index = bit2gf(xi_ind_bits,ind_bits_num) + 1;

                        x_rec_mat(xi_index,:) = xi_rec(1:end-ind_bits_num);

                        %x_rec_mat(xi_index,:) = xi_rec(Rseq(n_seg-k_seg+1:end));

                    end
                end

                %decoding in case of matched decoder detection
                if ind == 0;
                    %keep track of indexes that are already detected
                    detected_inds = zeros(1,seg_num);

                    for si = 1:seg_num,

                        metric_max = -1e10;

                        for sj = 1:seg_num,

                            if detected_inds(sj) == 0,

                                %calling the polar decoder
                                %[xi_temp, metric_j] = polar_dec(y_mat(si,deintliv_mat(sj,:)),k_seg,Rseq);

                                %out_rel_vec_intl = (1-p_bsc)*y_mat(si,deintliv_mat(sj,:))+(p_bsc)*(1-y_mat(si,deintliv_mat(sj,:)));

                                %[xi_temp, metric_j] = polar_dec_BSC(out_rel_vec_intl,k_seg,Rseq);

                                y_deliv = rem(y_mat(si,:)+coset_mat(sj,:),2);

                                out_rel_vec = (1-p_bsc)*y_deliv + (p_bsc)*(1-y_deliv);

                                [xi_temp, metric_j] = polar_dec_BSC(out_rel_vec,k_seg,Rseq);

                                %metric_j = metric_j

                                if metric_j > metric_max,

                                    metric_max = metric_j;

                                    xi_rec = xi_temp;

                                    xi_index = sj;

                                end


                            end

                        end

                        %xi_index = xi_index

                        detected_inds(xi_index) = 1;

                        x_rec_mat(xi_index,:) = xi_rec;

                    end
                end

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                %reshape polar decoders ordered output
                x_vec = reshape(transpose(x_rec_mat),1,prod(size(x_rec_mat)));

                %Remove the padding bits that were added at the encoder for segmentation
                x_vec = x_vec(1:n*m);

                %change the received bits from the polar decoder into symbols
                x_vec_symbols = bit2gf(x_vec,m);

                %Remove the padding bits that were added at the encoder for segmentation
                %x_vec_unpad = x_vec_symbols(1:n);

                %change the format of symbols into MATLAB's GF(2^m) symbols
                code_rec = gf(x_vec_symbols,m);

                %decode using RS decoder
                [rxcode,cnumerr] = rsdec(code_rec,n,k);

                rxint = zeros(1,k);

                for j=0:2^m-1,

                    rxint = rxint + j*(rxcode == j);

                end


                udec = gf2bit(rxint,m);

                %[u;udec]

                bit_errors = bit_errors + sum(u ~= udec);

                bx = sum(u ~= udec);

                if(bx > 0),
                    frame_errors = frame_errors+1;
                end

                symbol_errors = symbol_errors + sum(unb ~= rxint);

                %sb_ber = [sb, bit_errors/sb/lbits]

                %sb_fer = [sb, frame_errors]

                sb_ser = [sb, symbol_errors];

                sb_ind_k_p_fe = [sb, ind, k, p_bsc, frame_errors]

            end

            BER = bit_errors / lbits / sent_blocks;

            FER = frame_errors / sent_blocks;

            SER = symbol_errors / sent_blocks / length(rxint);

            if ind == 0,

                BER0(ik,ipb) = BER;
                FER0(ik,ipb) = FER;
                SER0(ik,ipb) = SER;

            end

            if ind == 1,

                BER1(ik,ipb) = BER;
                FER1(ik,ipb) = FER;
                SER1(ik,ipb) = SER;
                
            end


        end

    end

end

save BSC_results1.mat;


