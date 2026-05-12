
%IF Initial State
Y0 = inv(P0);
y0 = Y0*X0;
yhatu = zeros(nx,nt); yhatu(:,1) = y0;
Yu = zeros(nx,nx,nt); Yu(:,:,1) = Y0;

% define physical state for plotting
xinfo = zeros(nx, nt); xinfo(:,1) = X0;
Pinfo=zeros(nx,nx,nt); Pinfo(:,:,1) = P0;


% invert system matricies once for loop usage
Fi = inv(F); Qi = inv(Q); Ri = inv(R);

for k = 1:(nt-1)
    % predict 
    Mi = Fi'*Yu(:,:,k)*Fi;
    Om = Mi*G*inv(Qi + G'*Mi*G);
    Yp(:,:,k+1) = (eye(nx) - Om*G')*Mi;
    yhatp(:,K+1) = ((eye(nx) - Om*G')*Fi')*yhatu(:,k);

    % update
    yhatu(:,k+1) = yhatp(:,k+1) + H'*Ri*z(:,k+1);
    Yu(:,:,k+1) = Yp(:,:,k+1) + H'*Ri*H;

    % pull out state and covariance for plotting
    Pinfo(:,:,k+1) = inv(Yu(:,:,k+1));
    xinfo(:,k+1) = Pinfo(:,:,k+1)*yhatu(:,k+1);
end