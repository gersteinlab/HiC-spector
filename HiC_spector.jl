using HDF5;
using JLD;
using MAT;
using DataFrames;
using CurveFit;
using Interpolations;

function get_reproducibility(M1,M2,num_evec);

	if ~isequal(M1,M1')
		tmp1=M1-diagm(diag(M1));
		M1=tmp1+tmp1'+diagm(diag(M1));
	end
	if ~isequal(M2,M2')
		tmp2=M2-diagm(diag(M2));
		M2=tmp2+tmp2'+diagm(diag(M2));
	end

	k1=sum(sign(M1),2);
	d1=diag(M1);
	kd1=!((k1.==1).*(d1.>0))
	k2=sum(sign(M2),2);
	d2=diag(M2);
	kd2=!((k2.==1).*(d2.>0))
	iz=find((k1+k2.>0).*(kd1.>0).*(kd2.>0));

	M1b=M1[iz,iz];
	M2b=M2[iz,iz];

	i_nz1=find(sum(M1b,2).>0);
	i_nz2=find(sum(M2b,2).>0);

	i_z1=find(sum(M1b,2).==0);
	i_z2=find(sum(M2b,2).==0);

	Ln1_nz1=get_Laplacian(M1b);
	Ln2_nz2=get_Laplacian(M2b);

	(a1,b1)=eig(full(Ln1_nz1));
	(a2,b2)=eig(full(Ln2_nz2));

	ord1=sortperm(a1)[1:num_evec];
	b1=b1[:,ord1];
	ord2=sortperm(a2)[1:num_evec];
	b2=b2[:,ord2];

	b1_extend=zeros(size(M1b,1),num_evec);
	for i=1:num_evec
		b1_extend[i_nz1,i]=b1[:,i];
		x=b1_extend[:,i];
		x1=[x[2:end]' x[end]]';
		x2=[x[1] x[1:end-1]']';
		xx=(x+x1+x2)/3;
		b1_extend[i_z1,i]=xx[i_z1];
	end

	b2_extend=zeros(size(M2b,1),num_evec);
	for i=1:num_evec
		b2_extend[i_nz2,i]=b2[:,i];
		x=b2_extend[:,i];
		x1=[x[2:end]' x[end]]';
		x2=[x[1] x[1:end-1]']';
		xx=(x+x1+x2)/3;
		b2_extend[i_z2,i]=xx[i_z2];
	end

	evd=zeros(num_evec);
	for i=1:num_evec;
		evd[i]=evec_distance(b1_extend[:,i],b2_extend[:,i]);
	end

	Sd=sum(evd);

	evs=abs(sqrt(2)-Sd/num_evec)/sqrt(2);

	return evs,a1,a2;

end

function get_Laplacian(M);

	K=vec(sum(M,1));
	i_nz=find(K.>0);
	D_nz=spdiagm(K[i_nz]);
	D_isq=spdiagm(1./sqrt(K[i_nz]));

	Ln_nz=M[i_nz,i_nz]*D_isq;
	Ln_nz=I-D_isq*Ln_nz;
	n=size(M,1);

	Ln_nz=(Ln_nz+Ln_nz')/2;
	return Ln_nz;
end

function get_ipr(evec);

	ipr=1./sum(evec.^4,1);

end

function evec_distance(x,y);
	
	d1=sum((x-y).^2);
	d2=sum((x+y).^2);
	if d1<d2
		d=d1;
	else
		d=d2;
	end

	return sqrt(d);
end

#########################################################################################################################

#This code was adapted from the MATLAB code implemented in Knight and Ruiz, IMA Journal of Numerical Analysis (2012)
function knight_ruiz(M);
	M[isnan(M)]=0;
	L=size(M,1);
	iz=find(sum(M,2).>0);
	A=M[iz,iz];
	n=size(A,1);
	e = ones(n,1);
	res=[];
	delta = 0.1;
	x0 = e;
	tol = 1e-6;
	g=0.9; etamax = 0.1; # Parameters used in inner stopping criterion.

	eta = etamax;
	x = x0; rt = tol^2; v = x.*(A*x); rk = 1 - v;
	rho_km1=sum(rk.^2);
	rout = rho_km1; rold = rout;
	MVP = 0; # count matrix vector products.
	i = 0; # Outer iteration count.

	while rout > rt # Outer iteration
    	i = i + 1; k = 0; y = e;
    	innertol = maximum([eta^2*rout;rt]);
    	while rho_km1 > innertol #Inner iteration by CG
        	k = k + 1;
        	if k == 1
            	Z = rk./v; p=Z; rho_km1 = sum(rk.*Z);
        	else
            	beta=rho_km1/rho_km2;
            	p=Z + beta*p;
        	end
        	# Update search direction efficiently.
        	w = x.*(A*(x.*p)) + v.*p;
        	#w=squeeze(w,2);
        	alpha = rho_km1/sum(p.*w);
        	ap =squeeze(alpha*p,2);
        	# Test distance to boundary of cone.
        	ynew = y + ap;
        	if minimum(ynew) <= delta
            	if delta == 0
            		break
            	end
            	ind = find(ap .< 0);
            	gamma = minimum((delta - y[ind])./ap[ind]);
            	y = y + gamma*ap;
            	break
        	end
        	y = ynew;
        	rk = rk - alpha*w; rho_km2 = rho_km1; rho_km2=rho_km2[1];
        	Z = rk./v; rho_km1 = sum(rk.*Z);
    	end
    	x = x.*y; v = x.*(A*x);
    	rk = 1 - v; rho_km1 = sum(rk.*rk); rout = rho_km1;
    	MVP = MVP + k + 1;
    	# Update inner iteration stopping criterion.
    	rat = rout/rold; rold = rout; r_norm = sqrt(rout);
    	eta_o = eta; eta = g*rat;
    	if g*eta_o^2 > 0.1
        	eta = maximum([eta;g*eta_o^2]);
    	end
    	eta = maximum([minimum([eta;etamax]);0.5*tol/r_norm]);
    	#@sprintf("%3d %6d %.3e %.3e %.3e \n", i,k,r_norm,minimum(y),minimum(x));
        display(rout);
        #res=[res; r_norm];
	end
	#@printf("Matrix-vector products = %6d\n", MVP);
	x=squeeze(x,2);
	A2=A*diagm(x);
	A2=diagm(x)*A2;
	A_balance=extend_mat(A2,iz,L);
	A_balance=(A_balance+A_balance')/2;
	x_final=zeros(L);
	x_final[iz]=x;

	return x_final,A_balance;

end

function extend_mat(Z,iz,L);
    (u,v)=ind2sub(size(Z),find(Z.!=0));
    w=Z[find(Z)];
    #w=nonzeros(Z);
    u=iz[u];
    v=iz[v];
    Z_extend=sparse(u,v,w,L,L);
    Z_extend=full(Z_extend);
    return Z_extend;
end

#########################################################################################################################

function local_smoothing(x,y);

	span=0.01;
	v=sortperm(x);
	x=x[v];
	y=y[v];
	ux=unique(x);
	uy_smooth=zeros(size(ux));
	n=Int(floor(length(x)*span/2));

	mm=zeros(size(x));
	L=2*n+1;
	i=n+1;
	st=1;
	ed=i+n;
	mm[i]=mean(y[st:ed]);
	for i=n+2:length(y)-n;
		#display(i);
    	ed=ed+1;
    	mm[i]=mm[i-1]+y[ed]/L-y[st]/L;
	    st=st+1;
	end
	for i=1:n
    	mm[i]=mean(y[1:n+i]);
	end
	for i=1:n;
    	mm[end-n+i]=mean(y[end-n+1-n+i:end]);
	end

	for i=1:length(ux);
    	iz=find(x.==ux[i]);
    	uy_smooth[i]=mean(mm[iz]);
	end

	return ux,uy_smooth;

end

function get_expect_vs_d_single_chr_v0(W,chr2bins,bin_size);

	W=full(W);
	W[isnan(W)]=0;

	N=size(W,1);


	(u,v,w)=findnz(triu(W));
	d=float(v-u);
	d2=float(d);
	d2[d2.==0]=1/3;#this is the average distance for 2 points drawn from an uniform distribution between [0.1];
	d3=d2*bin_size;

	x=log10(d3);
	y=log10(w);

	xs,ys_smooth=local_smoothing(x,y);

	xs_all=collect(0:1.0:size(W,1)-1);xs_all[1]=1/3;
	xs_all=xs_all*bin_size;
	xs_all_aux=log10(xs_all);

	ys_all=zeros(size(xs_all));
	for k=1:length(xs_all_aux);
		ik=find(xs.==xs_all_aux[k]);
		if ~isempty(ik)
			ys_all[k]=ys_smooth[ik][1];
		end
	end

	A_x=find(ys_all.>0);
	knots=(A_x,);
	itp=interpolate(knots,ys_smooth, Gridded(Linear()));

	A_nz=find(ys_all.==0);
	for i=1:length(A_nz);
		ys_all[A_nz[i]]=itp[A_nz[i]];
	end

	expect=10.^ys_all;

	return xs_all, expect;

end

function get_expect_vs_d_WG_v0(contact,chr2bins,bin_size);

	all_d2=Float64[];
	all_w=Float64[];
	Ltmp=zeros(23);
	for chr_num=1:23

		#display(chr_num);
		W=extract_chr(contact,chr2bins,chr_num);
		W=full(W);
		W[isnan(W)]=0;

		N=size(W,1);

		(u,v,w)=findnz(triu(W));

		d=float(v-u);
		d2=float(d);
		d2[d2.==0]=1/3;

		all_d2=[all_d2;d2];
		all_w=[all_w;w];
		Ltmp[chr_num]=size(W,1);

	end

	all_d3=all_d2*bin_size;

	x=log10(all_d3);
	y=log10(all_w);

	xs,ys_smooth=local_smoothing(x,y);

	xs_all=collect(0:1.0:maximum(Ltmp)-1);xs_all[1]=1/3;
	xs_all=xs_all*bin_size;
	xs_all_aux=log10(xs_all);

	ys_all=zeros(size(xs_all));
	for k=1:length(xs_all_aux);
		ik=find(xs.==xs_all_aux[k]);
		if ~isempty(ik)
			ys_all[k]=ys_smooth[ik][1];
		end
	end

	A_x=find(ys_all.>0);
	knots=(A_x,);
	itp=interpolate(knots,ys_smooth, Gridded(Linear()));

	A_nz=find(ys_all.==0);
	for i=1:length(A_nz);
		ys_all[A_nz[i]]=itp[A_nz[i]];
	end

	expect=10.^ys_all;

	return xs_all, expect;

end


function extract_chr(A,chr2bins,chr_num);
	st=1+chr2bins[1,chr_num];
	ed=1+chr2bins[2,chr_num];
	A_chr=A[st:ed,st:ed];
	return A_chr;
end

function get_f_W(W,ys);

	N=size(W,1);
	W[isnan(W)]=0;
	dark_bins=find(sum(W,1).==0);
	num_dark=length(dark_bins);
	N_eff=N-num_dark;
	f_W=zeros(size(W));

	x=collect(1:N);

	for d=0:N-1
		f_W[1+d:N+1:end-d*N]=ys[d+1];
	end
	tmp=f_W-diagm(diag(f_W));
	f_W=f_W+tmp';
	#sum(f_W[1,:])=1 here..

	f_W[dark_bins,:]=0;
	f_W[:,dark_bins]=0;
	f_W=f_W/sum(f_W)*N_eff.^2;

	return f_W;

end

function get_compartment_A_B(W,f_W);

	iz=find(sum(W,2).>0);
	izz=find(sum(W,2).==0);
	Wn=W[iz,iz]./f_W[iz,iz];
	C=cor(Wn);
	(U,V)=eigs(C);
	i_max=indmax(U);
	ev=V[:,i_max:i_max+5];
	ev_whole=zeros(size(W,1),6);
	ev_whole[iz,:]=ev;
	ev_whole[izz,:]=NaN;

	(loc,span)=get_chunks_v2(sign(ev_whole[:,1]),1);#
	loc=round(loc);
	span=round(span);
	cpt=sign(ev_whole[loc,1]);

	return loc,span,ev_whole,cpt;

end

#id is the starting loc of a chunk, and d is the length it spans..
function get_chunks_v2(a,singleton=0);
	# adopt from a matlab code by Jiro Doke;
	 a                 = [NaN; a; NaN];
	 b                 = diff(a);
	 b1                = b;  # to be used in fullList (below)
	 ii                = trues(size(b));
	 ii[b.==0] = false;
	 b[ii]             = 1;
	 c                 = diff(b);
	 id                = find(c.==-1);

	 #Get single-element chunks also
	 if singleton.==1
	 	b1[id]          = 0;
	 	ii2             = find(b1[1:end-1]);
	 	d               = vcat(find(c.==1) - id + 1, ones(length(ii2)));
	 	id              = [id;ii2];
	 	v=sortperm(id);
	 	id=sort(id);
	 	#(id,tmp)        = sort(id);
	 	d               = d[v];
	 else
	 	d               = find(c.==1) - id + 1;
	 end

	 return id,d;
end

function report_compartments(hg19_info,bin2loc,loc,span,ev1,chr_num);
	
	iz=find(bin2loc[1,:].==chr_num-1);
	a=round(bin2loc[2:3,iz]');
	st=a[loc];
	ed=a[loc+span]-1;
	output=DataFrame();
	chr=cell(size(st));
	for i=1:length(chr);
		chr[i]=change_chr(hg19_info,chr_num);
	end
	output[:chr]=chr;
	output[:start]=st;
	output[:end]=ed;
	x=ev1[loc];
	x[isnan(x)]=0;
	output[:compartment]=sign(x);

	return output;
end

#########################################################################################################################

function generate_arbitrary_mapping_files(hg19_info,bin_size);

	num_of_chromosomes=size(hg19_info,1);
	chr2bins=zeros(2,num_of_chromosomes);
	chr_length=hg19_info[:length];
	chr_num_bins=round(Int64,floor(chr_length/bin_size))+1
	#chr_num_bins=int(floor(chr_length/bin_size))+1;
	chr2bins[2,:]=cumsum(chr_num_bins)'-1;
	chr2bins[1,1]=0;
	chr2bins[1,2:end]=chr2bins[2,1:end-1]+1;
	X=round(Int,chr2bins+1);
	bin2loc=zeros(3,X[2,end]);
	for c=1:size(hg19_info,1);
		bin2loc[1,X[1,c]:X[2,c]]=c-1;
		bin2loc[2,X[1,c]:X[2,c]]=round(Int,collect(1:bin_size:chr_length[c]))';
		bin2loc[3,X[1,c]:X[2,c]]=[round(Int,collect(bin_size:bin_size:chr_length[c]))' chr_length[c]];
	end
	return round(Int64,chr2bins),round(Int64,bin2loc);

end

function define_hg19_genome();

	hg19_info=DataFrame();
	hg19_info[:id]=1:25;
	hg19_info[:chr]=["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10",
	"chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21",
	"chr22","chrX","chrY","chrM"];
	hg19_info[:length]=[249250621,243199373,198022430,191154276,180915260,171115067,159138663,
	146364022,141213431,135534747,135006516,133851895,115169878,107349540,102531392,90354753,
	81195210,78077248,59128983,63025520,48129895,51304566,155270560,59373566,16571];

	return hg19_info;

end

function change_chr(hg19_info,chr)

	if typeof(chr)==Float64||typeof(chr)==Int64;
		chr2=hg19_info[:chr][hg19_info[:id].==chr][1];
	elseif typeof(chr)==ASCIIString||typeof(chr)==SubString{ASCIIString}||typeof(chr)==UTF8String
		chr2=hg19_infp[:id][hg19_info[:chr].==chr][1];
	end

	return chr2;

end


#input file is generated by HiC-Pro, with 3 columns: row_index, column_index, entry
function read_simple_contact_map(input_file,hg19_info,bin_size);

	chr_length=hg19_info[:length];
	X=readtable(input_file,separator='\t',header=false);
	chr2bins,bin2loc=generate_arbitrary_mapping_files(hg19_info,bin_size);
	N=size(bin2loc,1);
	M=sparse(X[:,1],X[:,2],X[:,3],N,N);
	if ~isequal(M,M');
		tmp=M-spdiagm(diag(M));
		M=M+tmp';
	end

	return M;
end

function define_ce10_genome();

	ce10_info=DataFrame();
	ce10_info[:id]=1:7;
	ce10_info[:chr]=["chrV","chrX","chrIV","chrII","chrI","chrIII","chrM"];
	ce10_info[:length]=[20924149,17718866,17493793,15279345,15072423,13783700,13794];
	return ce10_info;

end

function define_dm3_genome();

	dm3_info=DataFrame();
	dm3_info[:id]=1:15;
	dm3_info[:chr]=["2L","2LHet","2R","2RHet","3L","3LHet","3R","3RHet","4","U","Uextra","X","Xhet","YHet","dmel_mitochondrion_genome"];
	dm3_info[:length]=[23011544,368872,21146708,3288761,24543557,2555491,27905053,2517507,1351857,10049037,29004656,22422827,204112,347038,19517];

	return dm3_info;

end










