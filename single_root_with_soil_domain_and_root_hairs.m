%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Finite difference approach - single root with soil domain and root hairs
clc ; clear all;
load('Meunier_Conductances.mat') % parametrization: inverse modelling & root pressure probe       
%%
for i = 1:26
    pot(i) = {char(64+i)};
    pot(i+26)={char(96+i)};
end

psi_ss = linspace(-100,-10e3,10);   %soil water potential
Psi_collar = -15000;                %collar potential
rh_l = [0 0.04];                        %root hair length (Duddek et al. 2022)
rh_r = 17.88*1e-4;                  %root hair radius (Duddek et al. 2022)
r_z = 1;                             %rhi


%Sand
k0(1) = 0.047535;   %cm/s
h0(1)=  -3.402;       %this is alfa [cm^-1]
tau(1)= 3.6;
%Loam
k0(2)=  6.0e-4;       %cm/s
h0(2)=  -29.06;       %this is alfa [cm^-1]
tau(2)= 3;

ty = fieldnames(kr);
st = {'Sand' 'Loam'};
den = {'Hairless' 'Hairs'};

l = 0.1;    %root segment length


[Kx] = deal(cell(1, length(fieldnames(kr))));
[K_soil, K_r_com, psi_x, Jx, Jr,Kr] = deal(cell(length(fieldnames(kr)), 2, 2));
[KKsoil, Kcom] = deal(cell(4,2,2,length(psi_ss)));
[Jrall, Jxall,k_soil,k_hair,psi_xx] = deal(cell(length(fieldnames(kr)), 2, 2, length(psi_ss)));
[Sanddata, Loamdata] = deal(cell(1, 4));
[JrSanddata, JrLoamdata] = deal(cell(4, length(psi_ss)));

for i1 = 1:length(fieldnames(kr))
    L = length(kr.(ty{i1}))*0.1;
    N = round(L/l);

    [K_soil{i1,:,:}, K_r_com{i1,:,:}, psi_x{i1,:,:}, Jx{i1,:,:}, Jr{i1,:,:}] = deal(zeros(N, N));
    [Jrall{i1,:,:,:}, Jxall{i1,:,:,:}] = deal(zeros(N, N));
end

for z = 1:length(psi_ss)
    psi_s = psi_ss(z);
    for t = 1:length(fieldnames(kr))
        L = length(kr.(ty{t}))*0.1;
        N = round(L/l);
        RZ = zeros(N,N);
        for i = 1:N
            RZ(i,1:i)= 1;
        end
        RZ = fliplr(RZ);
        %Building the matrix (Meunier et al. 2017, Appendix A.)
        Firstline = [-1,zeros(1,2*N - 1)];
        left = diag(ones(1,N));
        right = left;
        
        seg_num = 2:N;
        prev = seg_num - 1; % vector of parent segments

        for i=1:(N-1)
           left(prev(i),i+1) = -1;
        end
        
        IM = [Firstline;left,right];
    
        Kx{t} = kx.(ty{t})(1:N)/l;

         for s = 1:2
            
            for hd = 1:2
                disp(['Cycle z: ',num2str(z),'; Cycle t: ',num2str(t),'; Cycle s: ',num2str(s),'; Cycle hd ',num2str(hd)])
                for mp = 1:N
                    for n = 1:N
                        Kr{t,s,hd}(mp,n) = 2.*pi.*r.(ty{t}).*l.*kr.(ty{t})(n);
                        K_soil{t,s,hd}(mp,n) = (2*pi*l)./log(r_z/(r.(ty{t})+(rh_l(hd)*RZ(mp,n))))*k0(s)*((psi_s/h0(s))^-tau(s));
                        K_r_com{t,s,hd}(mp,n) = (1./Kr{t,s,hd}(mp,n)+1./K_soil{t,s,hd}(mp,n)).^-1;
                    end
                    C = IM*diag([-Kx{t},-K_r_com{t,s,hd}(mp,1:N)]')*IM';
                    
                    c = C(2:(N+1),2:(N+1));
                    d = (-K_r_com{t,s,hd}(mp,1:N)*psi_s)' - [Kx{t}(1)*Psi_collar, zeros(1,N-1)]';
                    
                    psi_x{t,s,hd}(mp,1:N) = c\d  ;
                    Jx{t,s,hd}(mp,1:N) = cumsum(fliplr(((psi_s-psi_x{t,s,hd}(mp,1:N)').*K_r_com{t,s,hd}(mp,1:N)')')) ;
                    Jr{t,s,hd}(mp,1:N) = fliplr(((psi_s-psi_x{t,s,hd}(mp,1:N)').*K_r_com{t,s,hd}(mp,1:N)')');
                    Jrall{t,s,hd,z}(mp,1:N) = Jr{t,s,hd}(mp,1:N);
                    Jxall{t,s,hd,z}(mp,1:N) = Jx{t,s,hd}(mp,1:N);
                    KKsoil{t,s,hd,z} = K_soil{t,s,hd};
                    Kcom{t,s,hd,z} = K_r_com{t,s,hd}(mp,1:N);
                    psi_xx{t,s,hd,z}(mp,1:N) = psi_x{t,s,hd}(mp,1:N);
                end
            end
         end
    end

    for t = 1:4
        for hd = 2
        Sanddata{t,hd}(:,z) = (Jx{t,1,hd}(:,end)-Jx{t,1,1}(:,end))./Jx{t,1,1}(:,end);
        Loamdata{t,hd}(:,z) = (Jx{t,2,hd}(:,end)-Jx{t,2,1}(:,end))./Jx{t,2,1}(:,end);
        JrSanddata{t,hd,z}  = (Jr{t,1,hd}(:,:)-Jr{t,1,1}(:,:));
        JrLoamdata{t,hd,z}  = (Jr{t,2,hd}(:,:)-Jr{t,2,1}(:,:));
        end
    end
end
save('Finite_Diff_Meunier_MOD3_S1_Lowtest.mat','Sanddata','Loamdata','KKsoil','k_hair','k_soil','Kcom','ty','st','den','pot','psi_xx','Jx','Jr', 'JrSanddata','JrLoamdata','Jxall','Jrall','Kx','Kr')

