clear ,clc , close all;
t0 = clock;
% Parameter Setting
N=3;
M=224;
L=10000;

% Algorithm
load data30
X = reshape(X3D,L,M)'; % 3D to 2D

%Generate missing data(Revise Here)
Y = X;
pixelcut = [1:1702 , 1767:2500 , 2710:2722 , 2800:3500 , 3520 , 3545:3596 , 3600 , 3710 , 3800:10000]; 
bandcut = [1:220 , 223:224];                                                                           
Y(bandcut , pixelcut) =0;

% Subset 
Y_Omega = X;
Y_Omega(: , pixelcut) = [];
 
% Obtain A by HyperCSI
[~ , L] = size(Y_Omega);        
[A, timeHyperCSI,a] = HyperCSI(Y_Omega,N);

% Subset
A_Omega = A;    
A_Omega(bandcut,:) = [];   
Y_Omega2 = X;
Y_Omega2(bandcut , :) = [];

% PCA Dimension Reduction of Y_Omega2
[Y_re , L] = PCA(Y_Omega2,N);

% Obtain Hyperplane from A_Omega
[h_hat , b_hat , alpha_hat] = Hyperplane(A_Omega, N);

% Obtain S
S = ( h_hat*ones(1,L)- b_hat'*Y_re   ) ./ ( (  h_hat - sum( b_hat.*alpha_hat )' ) *ones(1,L) );
S(S<0) = 0;

% Obtain Y
Y_Ans = A*S;

% Performance Comparison
time = etime(clock,t0);
percentage = Missing_Percentage(Y)
Y1 = reshape(Y_Ans',100,100,224);
error = Error_Calculation(Y1,X3D)

% plot
%  x = a(1,:); y = a(2,:);
%  scatter(x,y,'filled');

figure;
map1_est= reshape(S(1,:),100,100);
subplot(2,3,1);
imshow(map1_est);title('map 1 est');

map2_est= reshape(S(2,:),100,100);
subplot(2,3,2);
imshow(map2_est);title('map 2 est');
 
map3_est= reshape(S(3,:),100,100);
subplot(2,3,3);
imshow(map3_est);title('map 3 est');
 
subplot(2,3,4);
plot(A(:,1)); title('est signature 1');
axis([1 224 0 1]);

subplot(2,3,5);
plot(A(:,2)); title('est signature 2');
axis([1 224 0 1]);

subplot(2,3,6);
plot(A(:,3)); title('est signature 3');
axis([1 224 0 1]);

%% subprogram 1 Hyperplane
function [h_hat, b_hat, Xd] = Hyperplane(X,N);
%------------------------ Step 1 ------------------------
[Xd,L] = PCA(X,N); % dimension reduced data (Xd is (N-1)-by-L)

%------------------------ Step 2 ------------------------
for i = 1:N
    bi_tilde(:,i) = compute_bi(Xd,i,N); % obtain bi_tilde
end

r = (1/2)*norm(Xd(:,1)-Xd(:,2));
for i = 1:N-1
    for j = i+1:N
        dist_ai_aj(i,j) = norm(Xd(:,i)-Xd(:,j));
        if (1/2)*dist_ai_aj(i,j) < r
            r = (1/2)*dist_ai_aj(i,j); % compute radius of hyperballs
        end
    end
end
Xd_divided_idx = zeros(L,1);
radius_square = r^2;
for k = 1:N
    [IDX_alpha_i_tilde]= find( sum(  (Xd- Xd(:,k)*ones(1,L) ).^2,1  )  < radius_square );
    Xd_divided_idx(IDX_alpha_i_tilde) = k ; % compute the hyperballs
end

%------------------------ Step 3 ------------------------
for i = 1:N
    Hi_idx = setdiff([1:N],[i]);
    for k = 1:1*(N-1)
        Ri_k = Xd(:,( Xd_divided_idx == Hi_idx(k) ));
        [val idx] = max(bi_tilde(:,i)'*Ri_k);
        pi_k(:,k) = Ri_k(:,idx); % find N-1 affinely independent points for each hyperplane
    end
    b_hat(:,i) = compute_bi([pi_k Xd(:,i)],N,N);
    h_hat(i,1) = max(b_hat(:,i)'*Xd);
end
end

%% subprogram 2 Normal Vector of Hyperplane
function [bi] = compute_bi(a0,i,N)
Hindx = setdiff([1:N],[i]);
A_Hindx = a0(:,Hindx);
A_tilde_i = A_Hindx(:,1:N-2)-A_Hindx(:,N-1)*ones(1,N-2);
bi = A_Hindx(:,N-1)-a0(:,i);
bi = (eye(N-1) - A_tilde_i*(pinv(A_tilde_i'*A_tilde_i))*A_tilde_i')*bi;
bi = bi/norm(bi);
return;
end

%% subprogram 3 PCA
function [Xd , L] = PCA(X,N)
[M L ] = size(X);
d = mean(X,2);
U = X-d*ones(1,L);
[eV D] = eig(U*U');
C = eV(:,M-N+2:end);
Xd = C'*(X-d*ones(1,L)); 
return;
end

%% subprogram 4 HyperCSI
function [A_est, time, a] = HyperCSI(X,N)
t0 = clock;
%%PCA Dimension Reduction
[M L ] = size(X);
row_mean_pca = mean(X,2);  %mean
U(1:224,:) = X(1:224,:)-row_mean_pca;

q1_idx=0;
q2_idx=0;
[V,D] = eig(U*U');
for i=1:224
   if D(i,i) > q1_idx
       q1_idx =i;
   end
end
D(q1_idx,q1_idx)=0;
for i=1:224
   if D(i,i) > q2_idx
       q2_idx = i;
   end
end
C(:,1)=V(:,q1_idx);
C(:,2)=V(:,q2_idx);
a = C'*U;
%=======================================
%%SPA Find the vertex
%Find frist vertex
tempnorm=0;
a_idx=0;
for i=1:L
    if norm(a(:,i)) > tempnorm
        tempnorm = norm(a(:,i));
        a_idx = i;
    end
end
a_spa(:,1) = a(:,a_idx);

%Find Second vertex
%generate projection matrix
PM = eye(2)-a_spa(:,1)*a_spa(:,1)' / (norm(a_spa(:,1))*norm(a_spa(:,1)));

tempnorm=0;
for j=1:L
    a_p(:,j) = PM*a(:,j);
    if norm(a_p(:,j)) > tempnorm
        tempnorm = norm(a_p(:,j));
        a_idx = j;
    end
end
a_spa(:,2) = a(:,a_idx);
a2_idx=a_idx;

%Find Third vertex
tempnorm=0;
for k=1:L
    a_p(:,k) = a_p(:,k)-a_p(:,a2_idx);
    if norm(a_p(:,k)) > tempnorm
        tempnorm = norm(a_p(:,k));
        a_idx = k;
    end
end
a_spa(:,3) = a(:,a_idx);
%=============================================
%Normal Vector Estimation
vd1 = a_spa(:,2)-a_spa(:,3);
vd2 = a_spa(:,3)-a_spa(:,1);
vd3 = a_spa(:,1)-a_spa(:,2);
b1_tilde = nv(vd1);
b2_tilde = nv(vd2);
b3_tilde = nv(vd3);
r1=norm(a_spa(:,2)-a_spa(:,3))*0.3;
r2=norm(a_spa(:,3)-a_spa(:,1))*0.3;
r3=norm(a_spa(:,1)-a_spa(:,2))*0.3;
max=0;max2=0;
%h1_hat of Normal Vector
for i=1:L
    v2 = a(:,i)-a_spa(:,2);
    v3 = a(:,i)-a_spa(:,3);
    if  norm(v2)<r2                        %testing point in given circle range
        if dot(a(:,i),b1_tilde)>max
            max = dot(a(:,i),b1_tilde);
            p1_idx1=i;
        end
    end
    if norm(v3)<r3                          %testing point in given circle range
        if dot(a(:,i),b1_tilde)>max2
            max2 = dot(a(:,i),b1_tilde);
            p2_idx1=i;
        end
    end
end
p1=a(:,p1_idx1)-a(:,p2_idx1);
b1_hat=nv(p1);

h1_hat=0;
for i=1:L
    m=b1_hat'*a(:,i);
    if m>h1_hat
        h1_hat = m;
    end
end
h1_hat=h1_hat/norm(b1_hat);
%h2_hat of Normal Vector
max=0;max2=0;
for i=1:L
    v3 = a(:,i)-a_spa(:,3);
    v1 = a(:,i)-a_spa(:,1);
    if norm(v3)<r3
        if dot(a(:,i),b2_tilde)>max
            max = dot(a(:,i),b2_tilde);
            p1_idx2=i;
        end
    end
    if norm(v1)<r1
        if dot(a(:,i),b2_tilde)>max2
            max2 = dot(a(:,i),b2_tilde);
            p2_idx2=i;
        end
    end
end
p2=a(:,p1_idx2)-a(:,p2_idx2);
b2_hat=nv(p2);
h2_hat=0;
for i=1:L
    m=b2_hat'*a(:,i);
    if m>h2_hat
        h2_hat = m;
    end
end
h2_hat=h2_hat/norm(b2_hat);
%h3_hat of Normal Vector
max=0;max2=0;
for i=1:L
    v1 = a(:,i)-a_spa(:,1);
    v2 = a(:,i)-a_spa(:,2);
    if  norm(v1)<r1
        if dot(a(:,i),b3_tilde)>max
            max = dot(a(:,i),b3_tilde);
            p1_idx3=i;
        end
    end
    if norm(v2)<r2
        if dot(a(:,i),b3_tilde)>max2
            max2 = dot(a(:,i),b3_tilde);
            p2_idx3=i;
        end
    end
end
p3=a(:,p1_idx3)-a(:,p2_idx3);
b3_hat=nv(p3);

h3_hat=0;
for i=1:L
    m=b3_hat'*a(:,i);
    if m>h3_hat
        h3_hat = m;
    end
end
h3_hat=h3_hat/norm(b3_hat);
%Find the Exact Vertex by Solving Simultaneous Equations
c1=h1_hat*norm(b1_hat);
c2=h2_hat*norm(b2_hat);
c3=h3_hat*norm(b3_hat);
A1(1,:)=b2_hat';
A1(2,:)=b3_hat';
A2(1,:)=b3_hat';
A2(2,:)=b1_hat';
A3(1,:)=b1_hat';
A3(2,:)=b2_hat';
B1=[c2,c3]';
B2=[c3,c1]';
B3=[c1,c2]';
P1=inv(A1)*B1;
P2=inv(A2)*B2;
P3=inv(A3)*B3;
%PCA inverse operation
A_est(:,1)=C*P1+row_mean_pca;
A_est(:,2)=C*P2+row_mean_pca;
A_est(:,3)=C*P3+row_mean_pca;
time = etime(clock,t0);
end

%%  subprogram 5 Normal Vector
function [new] = nv(old)
new(1,1) = -old(2,1);
new(2,1) = old(1,1);
return;
end