clc
clear all

filename = input('Enter the name of the text file: ','s');
fileID = fopen(filename,'r');
line = fgetl(fileID);
n = sscanf(line,'%f');
A = zeros(n,n);
for i=1:n
        line = fgetl(fileID);
        A(i,1:n) = sscanf(line,'%f');
end

line = fgetl(fileID);
max_iter = sscanf(line,'%f');

line = fgetl(fileID);
max_err = sscanf(line,'%f');

line = fgetl(fileID);
scaling_factor = sscanf(line,'%f');
    
disp('List of methods:');
    disp('1. Direct power method (to find the eigenvalue having the maximum magnitude and the corresponding eigenvector)  ');
    disp('2. Inverse power method (to find the eigenvalue having the minimum magnitude and the corresponding eigenvector)');
    disp('3. Shifted-power method (to find intermediate eigenvalues [based on Gershgorin disc] and corresponding eigenvectors) ');
    disp('4. QR method (to find all eigenvalues of a matrix; write your own program for GramSchmidt process)');

method_select = input('Select the method: ');

if(method_select == 1)
    x = ones(n,1);
    
    err = 100;
    iter = 0;
    eig = 0;
    while ((err > max_err) && (iter < max_iter))
        y = A*x
        s = max(abs(y(1:n,1)))
        x = y/s
        err = abs(100*(s - eig)/s);
        eig = s;
        iter = iter + 1;
    end
    x = x/norm(x);
    fileoID = fopen('Output1.txt','w');
    
    fprintf(fileoID,'Eigenvalue: ');
    fprintf(fileoID,'%.4f ',eig);
    fprintf(fileoID,'\r\n\n');
    fprintf(fileoID,'Eigenvector:\n');
    for i=1:n
        fprintf(fileoID,'%.4f  ',x(i,1));
        fprintf(fileoID,'\r\n');
    end
    fprintf(fileoID,'\r\n');
    fprintf(fileoID,'Iterations: ');
    fprintf(fileoID,'%d ',iter);
    fclose(fileID);
    fclose(fileoID);
end

if(method_select == 2)
    A = inv(A);
    x = ones(n,1);
    
    err = 100;
    iter = 0;
    eig = 0;
    while ((err >= max_err) && (iter < max_iter))
        y = A*x
        s = max(abs(y(1:n,1)))
        x = y/s
        err = abs(100*(1/s - eig)/1/s);
        eig = 1/s;
        iter = iter + 1;
    end
    x = x/norm(x);
    fileoID = fopen('Output2.txt','w');
    
    fprintf(fileoID,'Eigenvalue: ');
    fprintf(fileoID,'%.4f ',eig);
    fprintf(fileoID,'\r\n\n');
    fprintf(fileoID,'Eigenvector:\n');
    for i=1:n
        fprintf(fileoID,'%.4f  ',x(i,1));
        fprintf(fileoID,'\r\n');
    end
    fprintf(fileoID,'\r\n');
    fprintf(fileoID,'Iterations: ');
    fprintf(fileoID,'%d ',iter);
    fclose(fileID);
    fclose(fileoID);
end

if(method_select == 3)
    I = eye(n,n);
    sI = scaling_factor*I;
    B = A - sI;
    B = inv(B);
    x = ones(n,1);
    
    err = 100;
    iter = 0;
    eig = 0;
    while ((err > max_err) && (iter < max_iter))
        y = B*x
        s = max(abs(y(1:n,1)))
        x = y/s
        err = abs(100*((1/s - eig)/1/s));
        eig = 1/s;
        iter = iter + 1;
    end
    x = x/norm(x);
    fileoID = fopen('Output3.txt','w');
    
    fprintf(fileoID,'Eigenvalue: ');
    fprintf(fileoID,'%.4f ',eig+scaling_factor);
    fprintf(fileoID,'\r\n\n');
    fprintf(fileoID,'Eigenvector:\n');
    for i=1:n
        fprintf(fileoID,'%.4f  ',x(i,1));
        fprintf(fileoID,'\r\n');
    end
    fprintf(fileoID,'\r\n');
    fprintf(fileoID,'Iterations: ');
    fprintf(fileoID,'%d ',iter);
    fclose(fileID);
    fclose(fileoID);
    
end

if(method_select == 4)
    err = 100;
    c1=0;
    c2=0;
    for iter = 1:max_iter
        Q=[];
        a=[];
        t=[];
        for j=1:n
            a = A(:,j);
            for i=1:j-1
                a = a - dot(Q(:,i),A(:,j))*Q(:,i);
            end
            t = a/norm(a);
            Q = [Q, t];
        end

        r = Q' * A;
        A=r*Q;
        c2 = max(diag(A));
        if (iter>1)
            err = (100*(abs(c2-c1)/abs(c1)));
            if err < max_err
               break;
            end
        end
        c1=c2;
    end

    v = diag(A);
    fileoID = fopen('Output4.txt','w');
    fprintf(fileoID,'Eigenvalues:\n');
    for i=1:n
        fprintf(fileoID,'%.4f  ',v(i));
        fprintf(fileoID,'\r\n');
    end
    fprintf(fileoID,'\r\n');
    fprintf(fileoID,'Iterations: ');
    fprintf(fileoID,'%d ',iter);
    fclose(fileID);
    fclose(fileoID);
    disp(v);
end