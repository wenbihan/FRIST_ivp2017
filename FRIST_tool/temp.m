N = 10000;
n = 16;
sp = 0.02*n*N;
numIter = 1000;

W = dctmtx(sqrt(n));
W = kron(W, W);

A = randn(n, n);
[A, ~] = qr(A);

W0 = W;
W = A;
A = W0;
Atrans = A';


idx = randperm(n * N, sp);
X0 = zeros(n*N, 1);
X0(idx) = randn(sp, 1);
% X0(idx) = 1;
X0 = reshape(X0, [n, N]);

Y = Atrans * X0;

obj = zeros(numIter, 1);
relaChange = zeros(numIter, 1);

for iter = 1 : numIter
    
    W1 = W;
    X1 = W * Y;
    [s]=sort(abs(X1(:)),'descend');
%     X = X1.*(bsxfun(@ge, abs(X1(:)), s(sp)));
    X1(abs(X1) < s(sp)) = 0;
    X = reshape(X1, [n, N]);
    
    %Transform Update Step
    [U,~,V]=svd(Y*X');
    W=(V(:,1:n))*U';
    
    obj(iter) = norm(W * Y - X, 'fro')^2;
    relaChange(iter) = norm(W - W1, 'fro');
end

norm(W * Y - X, 'fro')
figure(1); plot(1:numIter, obj);
figure(3); plot(1:numIter, relaChange);


figure(2);imshow(abs(W * A'), [])