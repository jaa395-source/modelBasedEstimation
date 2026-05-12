%% MAE 6760 Model Based Estimation
% Cornell University
% M Campbell
%
% Homework #2
% Problem #5: Nonlinear Least Squares
%   for multi-sensor range localization problem
%
rng(10);

%% Problem set up with beacons and msmt generation
%Beacon locations
bA = [-10;100];
bB = [490;20];
bC = [500;40];
beacons=[bA bB bC];
%True location of the car
xtrue = [-5;2];

%Perfect measurements
RA = sqrt((bA(1) - xtrue(1))^2 + (bA(2) - xtrue(2))^2);
RB = sqrt((bB(1) - xtrue(1))^2 + (bB(2) - xtrue(2))^2);
RC = sqrt((bC(1) - xtrue(1))^2 + (bC(2) - xtrue(2))^2);

%a total of 10 measurements
n=10;


%% Part (a) - two ranging beacons (A,B,C)
%generate noisy measurements to beacons A,B,C
xhat = [0; 0];
Rpart_a = diag([10 10 10]);
v_a = sqrtm(Rpart_a)*randn(3,n);
z_a = repmat([RA;RB;RC],[1,n]) + v_a;

NLS_pass = 1; J=[]; Jold=1; iter=0;
while NLS_pass
    iter = iter + 1;
    xj = xhat(:,iter);
    Rhat_A = sqrt([xj(1) - bA(1)]^2 + [xj(2) - bA(2)]^2);
    Rhat_B = sqrt([xj(1) - bB(1)]^2 + [xj(2) - bB(2)]^2);
    Rhat_C = sqrt([xj(1) - bC(1)]^2 + [xj(2) - bC(2)]^2);

    H = [(xj(1) - bA(1))/Rhat_A (xj(2) - bA(2))/Rhat_A;
         (xj(1) - bB(1))/Rhat_B (xj(2) - bB(2))/Rhat_B;
         (xj(1) - bC(1))/Rhat_C (xj(2) - bC(2))/Rhat_C];

    M1=0; M2=0; J(iter) = 0;

    for k=1:n
        M1 = M1 + H'*inv(Rpart_a)*H;
        ek = [z_a(:,k)-[Rhat_A; Rhat_B; Rhat_C]];
        M2 = M2+H'*inv(Rpart_a)*ek;
        J(iter) - J(iter) + 0.5*ek'*inv(Rpart_a)*ek;
    end
    dxhat = inv(M1)*M2;
    xhat(:,iter+1) = xj + dxhat;
    if abs(Jold-J(iter)) < 1E-3
        NLS_pass = 0;
    else
        Jold = J(iter);
    end
end

 PxhatF = inv(M1);
 xhatF = xhat(:,end);

 disp("--- Data for Part A, Using All 3 Beacons ---");
 disp("Iterations: " + iter)
 disp("X_hat estimate: " + string(xhatF));
 disp("Error Covariance: " + string(PxhatF))



%% Part (b) - two ranging beacons (A,B)
%Get noisy measurements to beacons A,B from part (a)
%simply so that the same noise is used
xhat = [0; 0];
ii_beacons=[1 2];
Rpart_b1 = Rpart_a(ii_beacons,ii_beacons); 
z_b1 = z_a(ii_beacons,:);


NLS_pass = 1; J=[]; Jold=1; iter=0;
while NLS_pass
    iter = iter + 1;
    xj = xhat(:,iter);
    Rhat_A = sqrt([xj(1) - bA(1)]^2 + [xj(2) - bA(2)]^2);
    Rhat_B = sqrt([xj(1) - bB(1)]^2 + [xj(2) - bB(2)]^2);
    Rhat_C = sqrt([xj(1) - bC(1)]^2 + [xj(2) - bC(2)]^2);

    H = [(xj(1) - bA(1))/Rhat_A (xj(2) - bA(2))/Rhat_A;
         (xj(1) - bB(1))/Rhat_B (xj(2) - bB(2))/Rhat_B];

    M1=0; M2=0; J(iter) = 0;

    for k=1:n
        M1 = M1 + H'*inv(Rpart_b1)*H;
        ek = [z_b1(:,k)-[Rhat_A; Rhat_B]];
        M2 = M2+H'*inv(Rpart_b1)*ek;
        J(iter) - J(iter) + 0.5*ek'*inv(Rpart_b1)*ek;
    end
    dxhat = inv(M1)*M2;
    xhat(:,iter+1) = xj + dxhat;
    if abs(Jold-J(iter)) < 1E-3
        NLS_pass = 0;
    else
        Jold = J(iter);
    end
end

 PxhatF_b_1 = inv(M1);
 xhatF_b_1 = xhat(:,end);

 disp("--- Data for Part B, Using Beacons A and B ---");
 disp("Iterations: " + iter)
 disp("X_hat estimate: " + string(xhatF_b_1));
 disp("Error Covariance: " + string(PxhatF_b_1))

%% Part (b) - two ranging beacons (B,C)
%Get noisy measurements to beacons B,C from part (a)
%simply so that the same noise is used
ii_beacons=[2 3];
Rpart_b2 = Rpart_a(ii_beacons,ii_beacons); 
z_b2 = z_a(ii_beacons,:);
x0=[0;0];

NLS_pass = 1; J=[]; Jold=1; iter=0;
while NLS_pass
    iter = iter + 1;
    xj = xhat(:,iter);
    Rhat_A = sqrt([xj(1) - bA(1)]^2 + [xj(2) - bA(2)]^2);
    Rhat_B = sqrt([xj(1) - bB(1)]^2 + [xj(2) - bB(2)]^2);
    Rhat_C = sqrt([xj(1) - bC(1)]^2 + [xj(2) - bC(2)]^2);

    H = [(xj(1) - bB(1))/Rhat_B (xj(2) - bB(2))/Rhat_B;
         (xj(1) - bC(1))/Rhat_C (xj(2) - bC(2))/Rhat_C];

    M1=0; M2=0; J(iter) = 0;

    for k=1:n
        M1 = M1 + H'*inv(Rpart_b2)*H;
        ek = [z_b2(:,k)-[Rhat_B; Rhat_C]];
        M2 = M2+H'*inv(Rpart_b2)*ek;
        J(iter) - J(iter) + 0.5*ek'*inv(Rpart_b2)*ek;
    end
    dxhat = inv(M1)*M2;
    xhat(:,iter+1) = xj + dxhat;
    if abs(Jold-J(iter)) < 1E-3
        NLS_pass = 0;
    else
        Jold = J(iter);
    end
end

 PxhatF_b_2 = inv(M1);
 xhatF_b_2 = xhat(:,end);

 disp("--- Data for Part B, Using Beacons B and C ---");
 disp("Iterations: " + iter)
 disp("X_hat estimate: " + string(xhatF_b_2));
 disp("Error Covariance: " + string(PxhatF_b_2))

%% Part (c) - perfect linearization
xhat = [-5; 2];
ii_beacons=[1 2];
Rpart_b1 = Rpart_a(ii_beacons,ii_beacons); 
z_b1 = z_a(ii_beacons,:);


NLS_pass = 1; J=[]; Jold=1; iter=0;
while NLS_pass
    iter = iter + 1;
    xj = xhat(:,iter);
    Rhat_A = sqrt([xj(1) - bA(1)]^2 + [xj(2) - bA(2)]^2);
    Rhat_B = sqrt([xj(1) - bB(1)]^2 + [xj(2) - bB(2)]^2);
    Rhat_C = sqrt([xj(1) - bC(1)]^2 + [xj(2) - bC(2)]^2);

    H = [(xj(1) - bA(1))/Rhat_A (xj(2) - bA(2))/Rhat_A;
         (xj(1) - bB(1))/Rhat_B (xj(2) - bB(2))/Rhat_B];

    M1=0; M2=0; J(iter) = 0;

    for k=1:n
        M1 = M1 + H'*inv(Rpart_b1)*H;
        ek = [z_b1(:,k)-[Rhat_A; Rhat_B]];
        M2 = M2+H'*inv(Rpart_b1)*ek;
        J(iter) - J(iter) + 0.5*ek'*inv(Rpart_b1)*ek;
    end
    dxhat = inv(M1)*M2;
    xhat(:,iter+1) = xj + dxhat;
    if abs(Jold-J(iter)) < 1E-3
        NLS_pass = 0;
    else
        Jold = J(iter);
    end
end

 PxhatF_c_1 = inv(M1);
 xhatF_c_1 = xhat(:,end);

 disp("--- Data for Part C, Perfect Linerization, Using Beacons A and B ---");
 disp("Iterations: " + iter);
 disp("X_hat estimate: " + string(xhatF_c_1));
 disp("Error Covariance: " + string(PxhatF_c_1));

%% Part (c) perfect linearization - two ranging beacons (B,C)
%Get noisy measurements to beacons B,C from part (a)
%simply so that the same noise is used
ii_beacons=[2 3];
Rpart_b2 = Rpart_a(ii_beacons,ii_beacons); 
z_b2 = z_a(ii_beacons,:);
x0=[0;0];

NLS_pass = 1; J=[]; Jold=1; iter=0;
while NLS_pass
    iter = iter + 1;
    xj = xhat(:,iter);
    Rhat_A = sqrt([xj(1) - bA(1)]^2 + [xj(2) - bA(2)]^2);
    Rhat_B = sqrt([xj(1) - bB(1)]^2 + [xj(2) - bB(2)]^2);
    Rhat_C = sqrt([xj(1) - bC(1)]^2 + [xj(2) - bC(2)]^2);

    H = [(xj(1) - bB(1))/Rhat_B (xj(2) - bB(2))/Rhat_B;
         (xj(1) - bC(1))/Rhat_C (xj(2) - bC(2))/Rhat_C];

    M1=0; M2=0; J(iter) = 0;

    for k=1:n
        M1 = M1 + H'*inv(Rpart_b2)*H;
        ek = [z_b2(:,k)-[Rhat_B; Rhat_C]];
        M2 = M2+H'*inv(Rpart_b2)*ek;
        J(iter) - J(iter) + 0.5*ek'*inv(Rpart_b2)*ek;
    end
    dxhat = inv(M1)*M2;
    xhat(:,iter+1) = xj + dxhat;
    if abs(Jold-J(iter)) < 1E-3
        NLS_pass = 0;
    else
        Jold = J(iter);
    end
end

 PxhatF_c_2 = inv(M1);
 xhatF_c_2 = xhat(:,end);

 disp("--- Data for Part C, Perfect Linerization, Using Beacons B and C ---");
 disp("Iterations: " + iter);
 disp("X_hat estimate: " + string(xhatF_c_2));
 disp("Error Covariance: " + string(PxhatF_c_2));

%% Part (d) - correlated noise model
Rpart_d = [10 0 0; 0 10 9; 0 9 10];
v_d = sqrtm(Rpart_d)*randn(3,n);
z_d = repmat([RA;RB;RC],[1,n]) + v_d;

xhat = [0; 0];

NLS_pass = 1; J=[]; Jold=1; iter=0;
while NLS_pass
    iter = iter + 1;
    xj = xhat(:,iter);
    Rhat_A = sqrt([xj(1) - bA(1)]^2 + [xj(2) - bA(2)]^2);
    Rhat_B = sqrt([xj(1) - bB(1)]^2 + [xj(2) - bB(2)]^2);
    Rhat_C = sqrt([xj(1) - bC(1)]^2 + [xj(2) - bC(2)]^2);

    H = [(xj(1) - bA(1))/Rhat_A (xj(2) - bA(2))/Rhat_A;
         (xj(1) - bB(1))/Rhat_B (xj(2) - bB(2))/Rhat_B;
         (xj(1) - bC(1))/Rhat_C (xj(2) - bC(2))/Rhat_C];

    M1=0; M2=0; J(iter) = 0;

    for k=1:n
        M1 = M1 + H'*inv(Rpart_d)*H;
        ek = [z_d(:,k)-[Rhat_A; Rhat_B; Rhat_C]];
        M2 = M2+H'*inv(Rpart_d)*ek;
        J(iter) - J(iter) + 0.5*ek'*inv(Rpart_d)*ek;
    end
    dxhat = inv(M1)*M2;
    xhat(:,iter+1) = xj + dxhat;
    if abs(Jold-J(iter)) < 1E-3
        NLS_pass = 0;
    else
        Jold = J(iter);
    end
end

 PxhatF_d = inv(M1);
 xhatF_d = xhat(:,end);


 disp("--- Data for Part D, Correlated Covariance, Using All Three Beacons ---");
 disp("Iterations: " + iter);
 disp("X_hat estimate: " + string(xhatF_d));
 disp("Error Covariance: " + string(PxhatF_d));
