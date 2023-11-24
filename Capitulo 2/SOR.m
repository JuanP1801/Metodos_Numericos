function [s, n, R] = SOR(Tol, niter, w)

    o2 = [601, 151];
    o4 = [2571, 441];
    p1 = [0, 0];
    p2 = [-116, 690];
    p3 = [-139, 1407];
    delta2 = deg2rad(-15);
    delta3 = deg2rad(-30);

    x0 = [deg2rad(10), 800, deg2rad(20), deg2rad(30)]';

    r1 = p1 - o2;
    r2 = p2 - o2;
    r3 = p3 - o2;
    A = [1300*(pi/2) (pi/2)*deg2rad(200) 0 0;
         1300 deg2rad(200) 0 0;
         0 (pi/2)*(deg2rad(200)+delta2) 1300*(pi/2) 0;
         0 (pi/2)*(deg2rad(200)+delta3) 0 1300*(pi/2)];

    b = -[r1(1) r1(2) r2(1) r3(1)]';

    c = 0;
    error = Tol + 1;
    D = diag(diag(A));
    L = -tril(A, -1);
    U = -triu(A, +1);

    while error > Tol && c < niter
        T = inv(D - w*L) * ((1 - w)*D + w*U);
        C = w * inv(D - w*L) * b;
        x1 = T * x0 + C;
        E(c + 1) = norm(x1 - x0, 'inf');
        error = E(c + 1);
        x0 = x1;
        c = c + 1;
    end

    if error < Tol
        s = x0;
        n = c;
        R = max(abs(eig(T)));
        fprintf('La aproximación de la solución es:\n');
        disp(s);
        fprintf('Número de iteraciones realizadas: %d\n', n);
        fprintf('Radio espectral de T: %f\n', R);
        fprintf('Tolerancia alcanzada: %f\n', Tol);
    else
        s = x0;
        n = c;
        R = max(abs(eig(T)));
        detA = det(A);
        normA = norm(A, 'fro');  % Norma de Frobenius de A
        fprintf('La norma de la matriz A es: %f\n', normA);
        fprintf('El determinante de la matriz A es: %f\n', detA);
        
        fprintf('El radio espectral de T es: %f\n', R);
    end
end
