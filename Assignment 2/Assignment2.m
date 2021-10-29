clc
clear all

disp('List of methods:');
    disp('1. Gauss elimination (GE; without pivoting) ');
    disp('2. GE (with pivoting)');
    disp('3. GE (with scaling and pivoting)');
    disp('4. LU decomposition by using GE (without pivoting)');
    disp('5. LU decomposition by using GE (with pivoting)');
    disp('6. LU decomposition by using Crout method (without pivoting)');
    disp('7. Cholesky decomposition (for symmetric positive definite matrix))');

method_select = input('Select the method: ');

if(method_select == 1)
    filename = input('Enter the name of the text file: ','s');
    fileID = fopen(filename,'r');
    line = fgetl(fileID);
    n = sscanf(line,'%f');
    A = zeros(n,n+1);
    for i=1:n
            line = fgetl(fileID);
            A(i,1:n+1) = sscanf(line,'%f');
    end
    x = zeros(n,1);
    for i=1:n-1
        for j=i+1:n
            m = A(j,i)/A(i,i);
            A(j,:) = A(j,:) - m.*A(i,:);
        end
    end
    x(n,:) = A(n,n+1)/A(n,n);
    for i=n-1:-1:1
        summ = 0;
        for j=i+1:n
            summ = summ + A(i,j)*x(j,:);
        end
         x(i,:) = (A(i,n+1) - summ)/A(i,i);
    end
    
    fileoID = fopen('Output1.txt','w');
    fprintf(fileoID,'The Unknowns x:\n');
    for i=1:n
        fprintf(fileoID,'%.6f  ',x(i));
        fprintf(fileoID,'\r\n');
    end
    fprintf(fileoID,'\r\n');
    fprintf(fileoID,'Elements of U:\n');
    for i=1:n
        fprintf(fileoID,'%.6f  ',A(i,1:n));
        fprintf(fileoID,'\r\n');
    end
    fclose(fileID);
    fclose(fileoID);
end

if(method_select == 2)
    filename = input('Enter the name of the text file: ','s');
    fileID = fopen(filename,'r');
    line = fgetl(fileID);
    n = sscanf(line,'%f');
    A = zeros(n,n+1);
    for i=1:n
            line = fgetl(fileID);
            A(i,1:n+1) = sscanf(line,'%f');
    end
    
    x = zeros(n,1);   
    I = eye(n);
    for i=1:n-1
        [m,p]=max(abs(A(i:n,i)));
        c=A(i,:);
        A(i,:)=A(p+i-1,:);
        A(p+i-1,:)=c;
        
        d = I(i,:);
        I(i,:)=I(p+i-1,:);
        I(p+i-1,:)=d;
        for j=i+1:n
            m = A(j,i)/A(i,i);
            A(j,:) = A(j,:) - m.*A(i,:);
        end
    end
   
    x(n) = A(n,n+1)/A(n,n);
    for i=n-1:-1:1
        summ = 0;
        for j=i+1:n
            summ = summ + A(i,j)*x(j,:);
        end
        x(i,:) = (A(i,n+1) - summ)/A(i,i);
    end
    fileoID = fopen('Output2.txt','w');
    fprintf(fileoID,'The Unknowns x:\n');
    for i=1:n
        fprintf(fileoID,'%.6f  ',x(i));
        fprintf(fileoID,'\r\n');
    end
    fprintf(fileoID,'\r\n');
    fprintf(fileoID,'Permutation Matrix:\n');
    for i=1:n
        fprintf(fileoID,'%.6f  ',I(i,1:n));
        fprintf(fileoID,'\r\n');
    end
    fprintf(fileoID,'\r\n');
    fprintf(fileoID,'Elements of U:\n');
    for i=1:n
        fprintf(fileoID,'%.6f  ',A(i,1:n));
        fprintf(fileoID,'\r\n');
    end
    fclose(fileID);
    fclose(fileoID);
end

if(method_select == 3)
    filename = input('Enter the name of the text file: ','s');
    fileID = fopen(filename,'r');
    line = fgetl(fileID);
    n = sscanf(line,'%f');
    A = zeros(n,n+1);

    for i=1:n
            line = fgetl(fileID);
            A(i,:) = sscanf(line,'%f');
    end
    for i=1:n
        S(i,1) = max(A(i,:));
        A(i,:) = A(i,:)/S(i,1);
    end
    x = zeros(n,1);   
    I = eye(n);
    for i=1:n-1
        [m,p]=max(abs(A(i:n,i)));
        c=A(i,:);
        A(i,:)=A(p+i-1,:);
        A(p+i-1,:)=c;
        
        d = I(i,:);
        I(i,:)=I(p+i-1,:);
        I(p+i-1,:)=d;
        for j=i+1:n
            m = A(j,i)/A(i,i);
            A(j,:) = A(j,:) - m.*A(i,:);
        end
    end
    
    x(n) = A(n,n+1)/A(n,n);
    for i=n-1:-1:1
        summ = 0;
        for j=i+1:n
            summ = summ + A(i,j)*x(j,:);
        end
        x(i,:) = (A(i,n+1) - summ)/A(i,i);
    end
    fileoID = fopen('Output3.txt','w');
    fprintf(fileoID,'The Unknowns x:\n');
    for i=1:n
        fprintf(fileoID,'%.6f  ',x(i));
        fprintf(fileoID,'\r\n');
    end
    fprintf(fileoID,'\r\n');
    fprintf(fileoID,'Permutation Matrix:\n');
    for i=1:n
        fprintf(fileoID,'%.6f  ',I(i,1:n));
        fprintf(fileoID,'\r\n');
    end
    fprintf(fileoID,'\r\n');
    fprintf(fileoID,'Elements of U:\n');
    for i=1:n
        fprintf(fileoID,'%.6f  ',A(i,1:n));
        fprintf(fileoID,'\r\n');
    end
    fclose(fileID);
    fclose(fileoID);
end

if(method_select == 4)
    filename = input('Enter the name of the text file: ','s');
    fileID = fopen(filename,'r');
    line = fgetl(fileID);
    n = sscanf(line,'%f');
    A = zeros(n,n+1);
    
    for i=1:n
            line = fgetl(fileID);
            A(i,1:n+1) = sscanf(line,'%f');
    end
    L=zeros(n);
    U=zeros(n);
    for i = 1:n
       L(i,i) = 1;
    end
    for k=1:n
        for j=k:n
            s2=0;
            for m=1:k-1
                s2=s2+L(k,m)*U(m,j);
            end
            U(k,j)=(A(k,j)-s2);
        end    
        for i=k+1:n
            s1=0;
            for m=1:k-1
                s1=s1+L(i,m)*U(m,k);
            end
            L(i,k)=(A(i,k)-s1)/U(k,k);
        end
        
    end
    
    x = zeros(n,1);
    y = zeros(n,1);
    y(1) = A(1,n+1)/L(1,1);
    for i=2:n
        summ = 0;
        for j=1:i
            summ = summ + L(i,j)*y(j,:);
        end
        y(i,:) = (A(i,n+1) - summ)/L(i,i);
    end
    x(n) = y(n,:)/U(n,n);
    for i=n-1:-1:1
        summ = 0;
        for j=i+1:n
            summ = summ + U(i,j)*x(j,:);
        end
        x(i,:) = (y(i,:) - summ)/U(i,i);
    end
    
    fileoID = fopen('Output4.txt','w');
    fprintf(fileoID,'The Unknowns x:\n');
    for i=1:n
        fprintf(fileoID,'%.6f  ',x(i));
        fprintf(fileoID,'\r\n');
    end
    fprintf(fileoID,'\r\n');
    fprintf(fileoID,'L Matrix:\n');
    for i=1:n
        fprintf(fileoID,'%.6f  ',L(i,1:n));
        fprintf(fileoID,'\r\n');
    end
    fprintf(fileoID,'\r\n');
    fprintf(fileoID,'U Matrix:\n');
    for i=1:n
        fprintf(fileoID,'%.6f  ',U(i,1:n));
        fprintf(fileoID,'\r\n');
    end
    fclose(fileID);
    fclose(fileoID);
end

if(method_select == 5)
    filename = input('Enter the name of the text file: ','s');
    fileID = fopen(filename,'r');
    line = fgetl(fileID);
    n = sscanf(line,'%f');
    A = zeros(n,n+1);
    for i=1:n
            line = fgetl(fileID);
            A(i,1:n+1) = sscanf(line,'%f');
    end
    
    x = zeros(n,1);   
    I = eye(n);
    for i=1:n-1
        [m,p]=max(abs(A(i:n,i)));
        c=A(i,:);
        A(i,:)=A(p+i-1,:);
        A(p+i-1,:)=c;
        
        d = I(i,:);
        I(i,:)=I(p+i-1,:);
        I(p+i-1,:)=d;
        for j=i+1:n
            m = A(j,i)/A(i,i);
            A(j,:) = A(j,:) - m.*A(i,:);
        end
    end
    
    L=zeros(n);
    U=zeros(n);
    for i = 1:n
       L(i,i) = 1;
    end
    for k=1:n
        for j=k:n
            s2=0;
            for m=1:k-1
                s2=s2+L(k,m)*U(m,j);
            end
            U(k,j)=(A(k,j)-s2);
        end    
        for i=k+1:n
            s1=0;
            for m=1:k-1
                s1=s1+L(i,m)*U(m,k);
            end
            L(i,k)=(A(i,k)-s1)/U(k,k);
        end
    end
    
    x = zeros(n,1);
    y = zeros(n,1);
    y(1) = A(1,n+1)/L(1,1);
    for i=2:n
        summ = 0;
        for j=1:i
            summ = summ + L(i,j)*y(j,:);
        end
        y(i,:) = (A(i,n+1) - summ)/L(i,i);
    end
    x(n) = y(n,:)/U(n,n);
    for i=n-1:-1:1
        summ = 0;
        for j=i+1:n
            summ = summ + U(i,j)*x(j,:);
        end
        x(i,:) = (y(i,:) - summ)/U(i,i);
    end
    
    fileoID = fopen('Output5.txt','w');
    fprintf(fileoID,'The Unknowns x:\n');
    for i=1:n
        fprintf(fileoID,'%.6f  ',x(i));
        fprintf(fileoID,'\r\n');
    end
    fprintf(fileoID,'\r\n');
    fprintf(fileoID,'Permutation Matrix:\n');
    for i=1:n
        fprintf(fileoID,'%.6f  ',I(i,1:n));
        fprintf(fileoID,'\r\n');
    end
    fprintf(fileoID,'\r\n');
    fprintf(fileoID,'L Matrix:\n');
    for i=1:n
        fprintf(fileoID,'%.6f  ',L(i,1:n));
        fprintf(fileoID,'\r\n');
    end
    fprintf(fileoID,'\r\n');
    fprintf(fileoID,'U Matrix:\n');
    for i=1:n
        fprintf(fileoID,'%.6f  ',U(i,1:n));
        fprintf(fileoID,'\r\n');
    end
    fclose(fileID);
    fclose(fileoID);
end

if(method_select == 6)
    filename = input('Enter the name of the text file: ','s');
    fileID = fopen(filename,'r');
    line = fgetl(fileID);
    n = sscanf(line,'%f');
    A = zeros(n,n+1);
    for i=1:n
            line = fgetl(fileID);
            A(i,1:n+1) = sscanf(line,'%f');
    end
    
    L=zeros(n);
    U=zeros(n);
    for i = 1:n
       U(i,i) = 1;
    end
    for k=1:n
        for i=k:n
            s1=0;
            for m=1:k-1
                s1=s1+L(i,m)*U(m,k);
            end
            L(i,k)=A(i,k)-s1;
        end
        for j=k+1:n
            s2=0;
            for m=1:k-1
                s2=s2+L(k,m)*U(m,j);
            end
            U(k,j)=(A(k,j)-s2)/L(k,k);
        end    
    end
    
    x = zeros(n,1);
    y = zeros(n,1);
    y(1) = A(1,n+1)/L(1,1);
    for i=2:n
        summ = 0;
        for j=1:i
            summ = summ + L(i,j)*y(j,:);
        end
        y(i,:) = (A(i,n+1) - summ)/L(i,i);
    end
    x(n) = y(n,:)/U(n,n);
    for i=n-1:-1:1
        summ = 0;
        for j=i+1:n
            summ = summ + U(i,j)*x(j,:);
        end
        x(i,:) = (y(i,:) - summ)/U(i,i);
    end
    
    fileoID = fopen('Output6.txt','w');
    fprintf(fileoID,'The Unknowns x:\n');
    for i=1:n
        fprintf(fileoID,'%.6f  ',x(i));
        fprintf(fileoID,'\r\n');
    end
    fprintf(fileoID,'\r\n');
    fprintf(fileoID,'L Matrix:\n');
    for i=1:n
        fprintf(fileoID,'%.6f  ',L(i,1:n));
        fprintf(fileoID,'\r\n');
    end
    fprintf(fileoID,'\r\n');
    fprintf(fileoID,'U Matrix:\n');
    for i=1:n
        fprintf(fileoID,'%.6f  ',U(i,1:n));
        fprintf(fileoID,'\r\n');
    end
end

if(method_select == 7)
    filename = input('Enter the name of the text file: ','s');
    fileID = fopen(filename,'r');
    line = fgetl(fileID);
    n = sscanf(line,'%f');
    A = zeros(n,n+1);
    for i=1:n
            line = fgetl(fileID);
            A(i,1:n+1) = sscanf(line,'%f');
    end
    B=A(:,n+1);
    A=A(1:n,1:n);
    L = zeros(n);
    for j=1:n
        s=0;
        for k=1:j-1
            s=s+L(j,k)^2;
        end
        L(j,j)=sqrt(A(j,j)-s);
        for i=j+1:n
            s1=0;
            for k=1:j-1
                s1=s1+L(i,k)*L(j,k);
            end
            L(i,j)=(A(i,j)-s1)/L(j,j);
        end            
    end  
    U = transpose(L);
    x = zeros(n,1);
    y = zeros(n,1);
    y(1) = B(1,1)/L(1,1);
    for i=2:n
        summ = 0;
        for j=1:i
            summ = summ + L(i,j)*y(j,:);
        end
        y(i,:) = (B(i,:) - summ)/L(i,i);
    end
    x(n) = y(n,:)/U(n,n);
    for i=n-1:-1:1
        summ = 0;
        for j=i+1:n
            summ = summ + U(i,j)*x(j,:);
        end
        x(i,:) = (y(i,:) - summ)/U(i,i);
    end
    
    fileoID = fopen('Output7.txt','w');
    fprintf(fileoID,'The Unknowns x:\n');
    for i=1:n
        fprintf(fileoID,'%.6f  ',x(i));
        fprintf(fileoID,'\r\n');
    end
    fprintf(fileoID,'\r\n');
    fprintf(fileoID,'Cholesky factor:\n');
    for i=1:n
        fprintf(fileoID,'%.6f  ',L(i,1:n));
        fprintf(fileoID,'\r\n');
    end
end
