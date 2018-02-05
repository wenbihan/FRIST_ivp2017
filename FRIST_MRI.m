function [IOut, paramsout]= FRIST_MRI(I1,Q1,paramsin, FRISTparam)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Inputs: 1) I1 : Complex-valued reference image obtained from fully-sampled k-space data (If only undersampled measurements are available, please provide 
%                the zero-filling reconstruction as the reference. The sampling mask will then be applied in the k-space space of the zero-filled result.)
%        2) Q1 : Sampling Mask for k-space
%        3) paramsin: Structure that contains the input parameters of the simulation. The various fields are as follows -
%                   - lambda0: parameter to set the weight on the negative log-determinat+Frobenius norm terms in the problem formulation
%                   - nu: Weight on the data fidelity term in the problem formulation
%                   - C: maximum allowed norm for the image in the formulation. In this implementation, we assume that C is large so that the constraint is automatically 
%                        satisfied without explicitly enforcing it during iterations. The output parameter `nmexceed' can be used to monitor any violations of the constraint.
%                   - num: Number of iterations of the TLMRI algorithm
%                   - numiter: Number of iterations within transform learning (i.e., iterations of sparse coding and transform update, with fixed image)
%                   - n: Patch size, i.e., Total number of pixels in a square patch
%                   - N: Number of training signals used in the transform learning step of the algorithm.
%                   - W0: initial transform for the algorithm
%                   - r: Patch Overlap Stride (this implementation works with r=1)
%                   - s: This is a vector of the same length as the number of TLMRI algorithm iterations, and contains the respective 
%                        sparsity fractions (i.e., fraction of non-zeros in the sparse code matrix) to be used in the TLMRI algorithm iterations.
%                   - cini: If set to 1, the initial reconstruction in the TLMRI algorithm is set to the zero-filling reconstruction. 
%                           For all other values of `cini', an external input for the initial reconstruction (estimate) is used (see next parameter). 
%                   - initrecon: An initial image reconstruction  (estimate). Please ensure that this image has intensities (magnitudes) approximately in the range [0 1].
%                   - ct: If set to 1, the code additionally outputs various performance metrics computed over the algorithm iterations. Otherwise, set to 0.
%                   - co: If set to 1, the code additionally outputs various algorithm convergence metrics computed over the algorithm iterations. Otherwise, set to 0.
%                   - nl: If set to 1, it indicates that the input data is normalized (i.e., the peak intensity value in the reference image is 1). 
%                         For any other value of `nl', the code automatically applies a normalization before the algorithm begins.
%       4) FRISTparam : parameters for FRIST learning
%Outputs:  1) IOut: Reconstructed MR image.
%          2) paramsout - Structure containing various outputs, and convergence or performance metrics 
%                         for the TLMRI algorithm. Many of these are vectors (whose entries correspond to values at each iteration).
%                 - transform : Final learnt transform
%                 - PSNR0 : PSNR of initial reconstruction (output only when ct is set to 1)
%                 - PSNR : PSNR of the reconstruction at each iteration of the TLMRI algorithm (output only when ct is set to 1)
%                 - HFEN : HFEN of the reconstruction at each iteration of the TLMRI algorithm (output only when ct is set to 1)
%                 - runtime: total execution time for the algorithm (output only when ct is set to 1)
%                 - itererror : norm of the difference between the reconstructions at successive iterations (output only when co is set to 1)
%                 - obj: objective function values at each iteration of algorithm (output only when co is set to 1)
%                 - sp: sparsification error (computed over all patches) at each iteration of algorithm (output only when co is set to 1)
%                 - reg: value of the transform learning regularizer in the objective at each algorithm iteration (output only when co is set to 1)
%                 - dfit: value of the data fidelity component of the objective at each algorithm iteration (output only when co is set to 1)
%                 - nmexceed: Vector whose values are zero except for iterations where the norm constraint on the image is violated. 
%                             For such iterations, the actual norm of the image is stored. (output only when co is set to 1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Initializing algorithm parameters
[aa,bb]=size(I1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
paramsin.nu = 1e6/(aa*bb);
paramsin.N = aa * bb;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
La2=(paramsin.nu)*(aa*bb); 
num=paramsin.num; 
num2=paramsin.numiter; 
n=paramsin.n; 
N=paramsin.N;
C=paramsin.C;
et=paramsin.s;
l0=paramsin.lambda0;
wd=paramsin.r;
W=paramsin.W0;
ct=paramsin.ct;
nl=paramsin.nl;
cini=paramsin.cini;
co=paramsin.co;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(nl~=1)
sc = (max(max(abs(I1))));  %peak image intensity
I1=I1/sc; %normalization
end

I5 = fftshift(fft2(ifftshift(I1)));   %simulate k-space of input
index=find(Q1==1);  %Find sampled k-space locations from sampling mask Q1
I2=(double(I5)).*(Q1);    %Apply the mask in k-space

if(cini==1)
I11=fftshift(ifft2(ifftshift(I2)));  % Zero-filled reconstruction
else
I11=paramsin.initrecon;
end

I11p=I11;

%initializing performance and convergence metrics
if(ct==1)
ittime=zeros(1,num);highfritererror=zeros(1,num);PSNR1=zeros(1,num);
end
if(co==1)
obj=zeros(1,num); sp=zeros(1,num);reg=zeros(1,num);dfit=zeros(1,num);nmexceed=zeros(1,num);itererror=zeros(1,num);
end

%TLMRI iterations
for kp=1:num
    tic
    
    Iiter=I11;    
    %Creating image patches (including wrap around patches)
    Ib= [I11 I11(:,1:(sqrt(n)-1));I11(1:(sqrt(n)-1),:) I11(1:(sqrt(n)-1),1:(sqrt(n)-1))];
    [TE,idx] = my_im2col(Ib,[sqrt(n),sqrt(n)],wd);
    
    if(kp==1)
        N2=size(TE,2); %total number of overlapping image patches
        [rows,cols] = ind2sub(size(Ib)-sqrt(n)+1,idx);
    end
    %FRIST learning iterations
    FRISTparam.iter = num2;
    FRISTparam.l0 = l0;
    FRISTparam.frac = et(kp)*n*N2;

%     FRISTparam.mergeK = 25;

    [YH, W, outputParam]= FRIST_learning_merge_l0penalty(TE, W, FRISTparam);
    permutation = outputParam.permutation;
    IDX = outputParam.IDX;    
    %reconstruction
    X1 = sparseSail0p(W*YH, FRISTparam.frac);
    ZZ = W' * X1;
    ZZ = rotateBack_merge(ZZ, IDX, permutation);
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Image Update
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    
    %%%%%%%%%%%Code corresponding to Step 4.b.i in the Algorithm Pseudocode in the SIIMS paper referenced above%%%%%%%%%%%
    IMoutR=zeros(size(Ib));
    IMoutI=zeros(size(Ib));
    bbb=sqrt(n);
 
    for jj = 1:10000:N2
        jumpSize = min(jj+10000-1,N2);
        block=reshape(real(ZZ(:,jj:jumpSize)),bbb,bbb,jumpSize-jj+1);blockc=reshape(imag(ZZ(:,jj:jumpSize)),bbb,bbb,jumpSize-jj+1);
        for ii  = jj:jumpSize
            col = cols(ii); row = rows(ii);
            IMoutR(row:row+bbb-1,col:col+bbb-1)=IMoutR(row:row+bbb-1,col:col+bbb-1)+block(:,:,ii-jj+1);
            IMoutI(row:row+bbb-1,col:col+bbb-1)=IMoutI(row:row+bbb-1,col:col+bbb-1)+blockc(:,:,ii-jj+1);
        end;
    end

    IMout=IMoutR + (0+1i)*IMoutI;
    
    IMout2=zeros(aa,bb);
    IMout2(1:aa,1:bb)=IMout(1:aa,1:bb);
    IMout2(1:(sqrt(n)-1),:)= IMout2(1:(sqrt(n)-1),:)+ IMout(aa+1:size(IMout,1),1:bb);
    IMout2(:, 1:(sqrt(n)-1))=IMout2(:, 1:(sqrt(n)-1)) + IMout(1:aa,bb+1:size(IMout,2));
    IMout2(1:(sqrt(n)-1),1:(sqrt(n)-1))= IMout2(1:(sqrt(n)-1),1:(sqrt(n)-1))+ IMout(aa+1:size(IMout,1),bb+1:size(IMout,2));
    
    I2=fftshift(fft2(ifftshift(IMout2)));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
    %%%%%%%%%%%Code corresponding to Step 4.b.iii in the Algorithm Pseudocode in the SIIMS paper referenced above%%%%%%%%%%%
%     Lb3=fft2(IMout2Z); Lb=real(Lb3);
    Lb=n*ones(aa,bb);               % from SPIE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    Lb=fftshift(Lb);
    Lb(index)= Lb(index) + (La2);
    I2(index)= (I2(index) + (La2)*I5(index));
    
    I2= I2./Lb;   %Updated k-space
    I11=fftshift(ifft2(ifftshift(I2)));  %Updated image
    
    time=toc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Compute various performance or convergence metrics    
    if(ct==1)
    ittime(kp)=time;
    PSNR1(kp)=20*log10((sqrt(aa*bb))*1/norm(double(abs(I11))-double(abs(I1)),'fro'));
    highfritererror(kp)=norm(imfilter(abs(I11),fspecial('log',15,1.5)) - imfilter(abs(I1),fspecial('log',15,1.5)),'fro');
    end
    
    if(co==1)
        if(norm(abs(I11),'fro')>C)
        nmexceed(kp)=norm(abs(I11),'fro');
        end
        itererror(kp)= (norm(Iiter - I11,'fro'));
        Ib= [I11 I11(:,1:(sqrt(n)-1));I11(1:(sqrt(n)-1),:) I11(1:(sqrt(n)-1),1:(sqrt(n)-1))];
        [TE] = my_im2col(Ib,[sqrt(n),sqrt(n)],wd);
        N2=size(TE,2);
        obj(kp)= (norm((W*TE) - X1,'fro'))^2  + (l0*(N2)*(-log(abs(det(W))) + 0.5*((norm(W,'fro'))^2))) + (La2/(aa*bb))*((norm(I2(index) - I5(index),'fro'))^2);
        sp(kp)=(norm((W*TE) - X1,'fro'))^2;
        reg(kp)= (l0*(N2)*(-log(abs(det(W))) + 0.5*((norm(W,'fro'))^2)));
        dfit(kp)=(La2/(aa*bb))*((norm(I2(index) - I5(index),'fro'))^2);
    end   
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Outputs
IOut=I11;
paramsout.transform=W;
paramsout.permutation = permutation;
paramsout.IDX = IDX;
paramsout.numIter = kp;
if(ct==1)
    paramsout.PSNR0=20*log10((sqrt(aa*bb))*1/norm(double(abs(I11p))-double(abs(I1)),'fro'));
    paramsout.PSNR=PSNR1;
    paramsout.HFEN=highfritererror;
    paramsout.runtime=sum(ittime);
end
if(co==1)
paramsout.itererror=itererror;
paramsout.obj=obj;
paramsout.sp=sp;
paramsout.reg=reg;
paramsout.dfit=dfit;
paramsout.nmexceed=nmexceed;
end
end