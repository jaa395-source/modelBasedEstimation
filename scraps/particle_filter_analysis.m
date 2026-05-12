function [xhatu,Pu] = particle_filter_analysis(N_part,nx, nk, x0, P0)
X_part = zeros(nx, N_part, nk);
X_part(:,:,1) = mvnrnd(x0, P0, N_part)';
W0 = ones(1,N_part)/N_part;
W = zeros(nk,N_part); W(1,:)=W0;

for k = 1:nk-1
    % start with the prior particles
    X_prior = X_part(:,:,k);

    % Prediction Step
    W_part = Qsq*randn(nw, N_part);
    X_pred = predict_state_carpose(X_prior, W_part, dt);

    % Update Step
    z_current = z(:, k+1);
    Z_hat = H*X_pred;

    Inn = z_current*ones(1,N_part) - Z_hat;

    % Calculate likelihood of weighting of each particle
    for ip=1:N_part
        L(1,ip) = exp(-0.5*Inn(:,ip)'*inv(R)*Inn(:,ip));
    end

    % Update Weights
    Wk_unnorm = W(k,:).*L;
    Wk = Wk_unnorm/sum(Wk_unnorm);

    xhatu(:,k+1) = sum(Wk.*X_part(:,:,k+1), 2);
    Err = [X_part(:,:,k+1) - xhatu(:,k+1)*ones(1,N_part)];
    Pu(:,:,k+1) = Err*diag(Wk)*Err';
end

end