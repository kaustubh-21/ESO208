clc
clear all

filename = input('Enter the name of the text file: ','s');
fileID = fopen(filename,'r');
n = fscanf(fileID, '%f\n', 1);
for i = 1:n
     x(i) = fscanf(fileID, '%f ', 1);
     y(i) = fscanf(fileID, '%f\n', 1);
end
line = fgetl(fileID);
for i = 1:n-1
    X(i) = fscanf(fileID,'%f\n', 1);
end
line = fgetl(fileID);
s0 = fscanf(fileID,'%f', 1);
sn = fscanf(fileID,'%f', 1);

disp('List of methods:');
    disp('1. Linear Spline');
    disp('2. Quadratic Spline');
    disp('3. Natural Cubic Spline');
    disp('4. Not-a-Knot Cubic Spline');
    disp('5. Periodic Cubic Spline');
    disp('6. Clamped Cubic Spline');
    
method_select = input('Select the method: ');

if(method_select == 1) %Linear Spline
    grid on
    hold on
    title('Linear Spline Interpolation')
    plot(x,y,'o')
    A = zeros(n-1,1);
    B = zeros(n-1,1);
   
    for i = 1:n-1
       A(i) = (y(i+1) - y(i))/(x(i+1)-x(i));
       B(i) = y(i) - A(i)*x(i);
    end
    
    Y = zeros(n-1,1);
    for i=1:n-1
        for j=1:n-1
            if (X(i) > x(j) && X(i) < x(j+1))
                break
            end
        end
        Y(i)=B(j)+A(j)*X(i);
    end
    
    resolution = 100;
    for i = 1:n-1
        f = @(x) A(i)*x + B(i);
        xf = linspace(x(i), x(i+1), resolution);
        fig = plot(xf,f(xf),'b');
    end
    saveas(fig,'Linear Spline.png');
    hold off
    
    fileoID = fopen('Output Linear.txt','w');
    fprintf(fileoID,'Linear Spline: ');
    fprintf(fileoID,'\r\n\n');
    for i = 1:n-1
        fprintf(fileoID,'%.4f  ',X(i));
        fprintf(fileoID,'%.4f',Y(i));
        fprintf(fileoID,'\r\n');
    end
    fclose(fileID);
    fclose(fileoID);
end

if(method_select == 2) %Quadratic Spline
    u = zeros(1,n);
    for i = 2:n
        u(i) = 2*((y(i) - y(i-1))/(x(i) - x(i-1)))-u(i-1);
    end
    c = zeros(1,n-1);
    for i = 1:n-1
        c(i) = y(i)+ (u(i)*(x(i+1) - x(i)))/2;
    end 
    a0 = zeros(1,n-1);a1 = zeros(1,n-1);a2 = zeros(1,n-1);
    for i = 1:n-1
        l2 = 1;
        for j = 1:2
            a = poly(x(i));
            l2 = conv(l2,a);
        end
        term1 = ((u(i+1)/(x(i+1) - x(i)))*l2)/2;
        l1 = 1;
        for j = 1:2
            a_1 = poly(x(i+1));
            l1 = conv(l1,a_1);
        end
        term2 = ((u(i)/(x(i+1) - x(i)))*l1)/2;
        p = term1 - term2;
        p(3) = p(3) + c(i);
        a0(i) = p(3);a1(i) = p(2);a2(i) = p(1);
    end
    
    fid = fopen('Output Quadratic.txt','w');
    l2 = n;
    for k = 1:n-1
        fprintf(fid,'%f %f %f in [%f,%f]\n',a2(k),a1(k),a0(k),x(k),x(k+1));
    end
    x1 = x(1):0.01:x(2);
    y1 = polyval([a2(1),a1(1),a0(1)],x1);
    plot (x1,y1,'b');
    hold on
    plot (x,y,'o');
    for i = 2:n-1
        x1 = x(i):0.01:x(i+1);
        y1 = polyval([a2(i),a1(i),a0(i)],x1);
        plot(x1,y1,'b');
    end
  
    hold off
end

if(method_select == 3) %Natural Cubic Spline
    
    h = zeros(n-1,1);
    for j = 1:n-1
        h(j) = x(j+1) - x(j);
    end

    g = zeros(n-1,1);
    for j = 1:n-1
        g(j) = (y(j+1) - y(j))/h(i);
    end
    
    s = zeros(n,1);
    a = zeros(n,n);
    q = zeros(n-2,1);
    for i = 1:n-2
       a(i+1,i) = h(i);
       a(i+1,i+1) = 2*(h(i)+h(i+1));
       a(i+1,i+2) = h(i+1);
    end
    for i = 1:n-2
        q(i,1) = 6*(g(i+1) - g(i));
    end
    p1 = a(2:n-1,2:n-1);
    r = inv(p1)*q;
    for i = 2:n-1
        s(i)= r(i-1);
    end
    
    A = zeros(n-1,1);
    B = zeros(n-1,1);
    C = zeros(n-1,1);
    D = zeros(n-1,1);
    for i = 1:n-1
       A(i) = s(i+1)/(6*h(i));
       B(i) = s(i)/(6*h(i));
       C(i) = y(i+1)/h(i) - A(i)*h(i)*h(i);
       D(i) = y(i)/h(i) - B(i)*h(i)*h(i);  
    end
    
    Y = zeros(n-1,1);
    for i=1:n-1
        for j=1:n-1
            if (X(i) > x(j) && X(i) < x(j+1))
                break
            end
        end
        Y(i)=A(i)*(X(i)-x(j))^3 - B(i)*(X(i)-x(j+1))^3 + C(i)*(X(i)-x(j)) - D(i)*(X(i)-x(j+1));
    end
    
    grid on
    hold on
    title('Natural Cubic Spline');
    plot(x,y,'o');
    for i = 1:n-1
       resolution=100;
        o=x(i);
        u=x(i+1);
        r=A(i);
        s=B(i);
        c=C(i);
        d=D(i);
        f = @(x) r*(x-o).*(x-o).*(x-o)-s*(x-u).*(x-u).*(x-u)+c*(x-o)-d*(x-u);
        xf = linspace(o, u, resolution);
        fig = plot(xf,f(xf),'b');
    end
    saveas(fig,'Natural Cubic Spline.png');
    hold off
    
    fileoID = fopen('Output Natural Spline.txt','w');
    fprintf(fileoID,'Natural Spline: ');
    fprintf(fileoID,'\r\n\n');
    for i = 1:n-1
        fprintf(fileoID,'%.4f  ',X(i));
        fprintf(fileoID,'%.4f',Y(i));
        fprintf(fileoID,'\r\n');
    end
    fclose(fileID);
    fclose(fileoID);
end

if(method_select == 4) %Not-a-Knot Cubic Spline
    
    h = zeros(n-1,1);
    for j = 1:n-1
        h(j) = x(j+1) - x(j);
    end

    a = zeros(n,n);
    a(1,1) = h(2,1);
    a(1,2) = -(h(1,1)+h(2,1));
    a(1,3) = h(1,1);
    a(n,n-2) = h(n-1,1);
    a(n,n-1) = -(h(n-1,1)+h(n-2,1));
    a(n,n) = h(n-2,1);
    for i = 2:n-1
        a(i,i-1) = h(i-1,1);
        a(i,i) = 2*(h(i-1,1)+h(i-1,1));
        a(i,i+1) = h(i-1,1);
    end
    
    b = zeros(n,1);
    g = zeros(n,1);
    for i=2:n
        g(i) = (y(i)-y(i-1))/h(i-1);
    end
    for i = 2:n-1
        b(i,1) = 6*(g(i+1,1)-g(i,1));
    end
    
    s = inv(a)*b;
    A = zeros(n-1,1);
    B = zeros(n-1,1);
    C = zeros(n-1,1);
    D = zeros(n-1,1);
    for i=1:n-1
        A(i,1)=s(i+1)/(6*h(i));
        B(i,1)=s(i)/(6*h(i));
        C(i,1)=(y(i+1)/h(i))-(s(i+1)*h(i)/6);
        D(i,1)=(y(i)/h(i))-(s(i)*h(i)/6);
    end
    
    
    Y = zeros(n-1,1);
    for i=1:n-1
        for j=1:n-1
            if (X(i) > x(j) && X(i) < x(j+1))
                break
            end
        end
        Y(i)=A(i)*(X(i)-x(j))^3 - B(i)*(X(i)-x(j+1))^3 + C(i)*(X(i)-x(j)) - D(i)*(X(i)-x(j+1));
    end
    
    resolution = 100;
    grid on
    hold on
    title('Not-a-Knot Cubic Spline');
    plot(x,y,'o');
    for i = 1:n-1
       resolution=100;
        o=x(i);
        u=x(i+1);
        r=A(i);
        s=B(i);
        c=C(i);
        d=D(i);
        f = @(x) r*(x-o).*(x-o).*(x-o)-s*(x-u).*(x-u).*(x-u)+c*(x-o)-d*(x-u);
        xf = linspace(o, u, resolution);
        fig = plot(xf,f(xf),'b');
    end
    saveas(fig,'Not-a-Knot Spline.png');
    hold off
       
    fileoID = fopen('Output Not-a-Knot Spline.txt','w');
    fprintf(fileoID,'Not-a-Knot Spline: ');
    fprintf(fileoID,'\r\n\n');
    for i = 1:n-1
        fprintf(fileoID,'%.4f  ',X(i));
        fprintf(fileoID,'%.4f',Y(i));
        fprintf(fileoID,'\r\n');
    end
    fclose(fileID);
    fclose(fileoID);
end

if(method_select == 5) %Preiodic Cubic Spline
   
    h=zeros(n-1);
    for i=1:n-1
        h(i)=x(i+1)-x(i);
    end

    a=zeros(n,n);
    a(1,1)=2*h(1);
    a(1,2)=h(1);
    a(1,n-1)=h(n-1);
    a(1,n)=2*h(n-1);
    a(n,1)=1;
    a(n,n)=-1;
    for i=2:n-1
        a(i,i-1)=h(i-1);
        a(i,i)=2*(h(i-1)+h(i));
        a(i,i+1)=h(i);
    end
    b=zeros(n);
    g=zeros(n);
    for i=2:n
        g(i) = (y(i)-y(i-1))/h(i-1);
    end
    for i=2:n-1
        b(i) = 6*(g(i+1)-g(i));
    end
    b(1)=-6*((y(n)-y(n-1))/h(n-1)) + (6*(y(2) - (y(1)))/h(1));
    s=inv(a)*b;
    A=zeros(n-1);
    B=zeros(n-1);
    C=zeros(n-1);
    D=zeros(n-1);
    for i=1:n-1
        A(i)=s(i+1)/(6*h(i));
        B(i)=s(i)/(6*h(i));
        C(i)=(y(i+1)/h(i))-(s(i+1)*h(i)/6);
        D(i)=(y(i)/h(i))-(s(i)*h(i)/6);
    end
    
    Y = zeros(n-1,1);
    for i=1:n-1
        for j=1:n-1
            if (X(i) > x(j) && X(i) < x(j+1))
                break
            end
        end
        Y(i)=A(i)*(X(i)-x(j))^3 - B(i)*(X(i)-x(j+1))^3 + C(i)*(X(i)-x(j)) - D(i)*(X(i)-x(j+1));
    end
    
    grid on
    hold on
    title('Periodic Cubic Spline');
    plot(x,y,'o');
    for i = 1:n-1
       resolution=100;
        o=x(i);
        u=x(i+1);
        r=A(i);
        s=B(i);
        c=C(i);
        d=D(i);
        f = @(x) r*(x-o).*(x-o).*(x-o)-s*(x-u).*(x-u).*(x-u)+c*(x-o)-d*(x-u);
        xf = linspace(o, u, resolution);
        fig = plot(xf,f(xf),'b');
    end
    saveas(fig,'Periodic Cubic Spline.png');
    hold off
    
    fileoID = fopen('Output Periodic.txt','w');
    fprintf(fileoID,'Periodic Spline: ');
    fprintf(fileoID,'\r\n\n');
    for i = 1:n-1
        fprintf(fileoID,'%.4f  ',X(i));
        fprintf(fileoID,'%.4f',Y(i));
        fprintf(fileoID,'\r\n');
    end
    fclose(fileID);
    fclose(fileoID);
end

if(method_select == 6) %Clamped Cubic Spline
    
    h=zeros(n-1,1);
    for i=1:n-1
        h(i)=x(i+1)-x(i);
    end
    a=zeros(n,n);
    a(1,1)=2*h(1);
    a(1,2)=h(1);
    a(n,n)=2*h(n-1);
    a(n,n-1)=h(n-1);
    for i=2:n-1
        a(i,i-1)=h(i-1);
        a(i,i)=2*(h(i-1)+h(i));
        a(i,i+1)=h(i);
    end
    b=zeros(n);
    g=zeros(n);
    for i=2:n
        g(i) = (y(i)-y(i-1))/h(i-1);
    end
    for i=2:n-1
        b(i) = 6*(g(i+1)-g(i));
    end
    b(1)=6*(((y(2)-y(1))/h(1))-s0);
    b(n)=6*(((y(n-1)-y(n))/h(n-1)) + sn);
    s=inv(a)*b;
    A=zeros(n-1);
    B=zeros(n-1);
    C=zeros(n-1);
    D=zeros(n-1);
    for i=1:n-1
        A(i)=s(i+1)/(6*h(i));
        B(i)=s(i)/(6*h(i));
        C(i)=(y(i+1)/h(i))-(s(i+1)*h(i)/6);
        D(i)=(y(i)/h(i))-(s(i)*h(i)/6);
    end
    
    Y = zeros(n-1,1);
    for i=1:n-1
        for j=1:n-1
            if (X(i) > x(j) && X(i) < x(j+1))
                break
            end
        end
        Y(i)=A(i)*(X(i)-x(j))^3 - B(i)*(X(i)-x(j+1))^3 + C(i)*(X(i)-x(j)) - D(i)*(X(i)-x(j+1));
    end
    
    grid on
    hold on
    title('Clamped Cubic Spline');
    plot(x,y,'o');
    for i = 1:n-1
       resolution=100;
        o=x(i);
        u=x(i+1);
        r=A(i);
        s=B(i);
        c=C(i);
        d=D(i);
        f = @(x) r*(x-o).*(x-o).*(x-o)-s*(x-u).*(x-u).*(x-u)+c*(x-o)-d*(x-u);
        xf = linspace(o, u, resolution);
        fig = plot(xf,f(xf),'b');
    end
    saveas(fig,'Clamped Cubic Spline.png');
    hold off
    
    fileoID = fopen('Output Clamped Spline.txt','w');
    fprintf(fileoID,'Clamped Cubic Spline: ');
    fprintf(fileoID,'\r\n\n');
    for i = 1:n-1
        fprintf(fileoID,'%.4f  ',X(i));
        fprintf(fileoID,'%.4f',Y(i));
        fprintf(fileoID,'\r\n');
    end
    fclose(fileID);
    fclose(fileoID);
end
    