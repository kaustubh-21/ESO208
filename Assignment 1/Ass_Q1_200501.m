clc
clear all

syms x

func = input('Give an equation in x: ','s')  ;
func = str2func(['@(x)',func]) ;

figure(1)
fplot(func)
title('f(x) vs x')
xlabel('x');
ylabel('f(x)');

disp('List of methods:');
    disp('1. Bisection Method');
    disp('2. False-Position');
    disp('3. Modified Flase-Position');
    disp('4. Newton Raphson');
    disp('5. Secant');
    
method_select = input('Select the method for finding Root: ') ;


if(method_select == 1)  %Bisection method
    x0 = input('Enter x0: ') ;
    x1 = input('Enter x1: ') ;
    max_err = input('Enter maximum Error: ') ;
    max_iter = input('Enter maximum Iteration allowed: ') ;
    if(func(x0)*func(x1)>0)
        disp('No root\n') ;
    else
        x_mid = (x0+x1)/2 ;
        error = 100 ;
        iter = 0 ;
        while ((error > max_err) && (iter < max_iter))
           if(func(x_mid) == 0)
               break;
           end
           if(func(x0)*func(x_mid) < 0)
               x1 = x_mid ;
           else
               x0 = x_mid ;
           end
           x_mid_new = (x0+x1)/2 ;
           error = abs(100*(x_mid_new-x_mid)/x_mid_new) ;
           err_vs_iter(iter + 1) = error ;
           x_mid = x_mid_new ;
           iter = iter + 1;
        end
        disp('Root is: ') ;
        disp(x_mid) ;
    end
    figure(2)
    plot(err_vs_iter)
    title('Bisection')
    xlabel('Iteration');
    ylabel('Error');
end

if(method_select == 2)   %False Position
    xl = input('Enter x0: ') ;
    xu = input('Enter x1: ') ;
    max_err = input('Enter maximum Error: ') ;
    max_iter = input('Enter maximum Iteration allowed: ') ;
    error = 100 ;
    iter = 0 ;
    if(func(xl)*func(xu) > 0)
        disp('No root') ;
    else
        xk = xl - (xl - xu)*func(xl)/(func(xl) - func(xu));
        
        while ((error > max_err) && (iter < max_iter))
            
            if(func(xl)*func(xu) < 0)
                xu = xk ;
            else
                xl = xk ;
            end
            x_k1 = xl - (xl - xu)*func(xl)/(func(xl) - func(xu));
            error = abs(100*(x_k1-xk)/x_k1) ;
            err_vs_iter(iter + 1) = error ;
            xk = x_k1 ;
            iter = iter + 1;
        end
        disp('Root is: ') ;
        disp(xk) ;
    end
    figure(2)
    plot(err_vs_iter)
    title('False Position')
    xlabel('Iteration');
    ylabel('Error');
end

if(method_select == 3)   %Modified False Position
    xl = input('Enter x0: ') ;
    xu = input('Enter x1: ') ;
    max_err = input('Enter maximum Error: ') ;
    max_iter = input('Enter maximum Iteration allowed: ') ;
    error = 100 ;
    iter = 0 ;
    f_low = func(xl) ;
    f_up = func(xu) ;
    if(f_low*f_up > 0)
        disp('No root') ;
    else
        xk = xl - (xl - xu)*func(xl)/(func(xl) - func(xu));
        while ((error >= max_err) && (iter < max_iter))
            x_k1 = xl - (xl - xu)*func(xl)/(func(xl) - func(xu)/2);
            if(func(xl)*func(xu) < 0)
                xu = x_k1 ;
            else
                xl = x_k1 ;

            end
            error = abs(100*(x_k1-xk)/x_k1) ;
            err_vs_iter(iter + 1) = error ;
            if(func(x_k1) == 0)
                break;
            end
            xk = x_k1 ;
            iter = iter + 1;
        end
        disp('Root is: ') ;
        disp(xk) ;
    end
    figure(2)
    plot(err_vs_iter)
    title('Modified False Position')
    xlabel('Iteration');
    ylabel('Error');
end

if(method_select == 4)   %Newton Rapson
    x0 = input('Enter x0: ') ;
    max_err = input('Enter maximum error allowed: ') ;
    max_iter = input('Enter Maximum Iterations allowed: ') ;
    func_der(x) = diff(func,x) ;
    disp(func_der);
    error = 100 ;
    iter = 0 ;
    
    while ((error > max_err) && (iter < max_iter))
        x1 = x0 - func(x0)/func_der(x0);
        error = abs(100*(x1-x0)/x1) ;
        err_vs_iter(iter + 1) = error ;
        iter = iter + 1 ;
        x0 = x1;
        double(x0);
    end
    disp('Root is: ') ;
    disp(vpa(x0,4));
    
    figure(2)
    plot(err_vs_iter)
    title('Newton Rapson')
    xlabel('Iteration');
    ylabel('Error');
end 

if(method_select == 5)   %Secant
    x0 = input('Enter x0: ') ;
    x1 = input('Enter x1: ') ;
    max_err = input('Enter maximum Error: ') ;
    max_iter = input('Enter maximum Iteration allowed: ') ;
    error = 100 ;
    iter = 0 ;
    
    while ((error > max_err) && (iter < max_iter))
        x2 = x1 - func(x1)*(x1-x0)/(func(x1)-func(x0)) ;
        error = abs(100*(x2-x1)/x2) ;
        err_vs_iter(iter + 1) = error ;
        x0 = x1 ;
        x1 = x2 ;
        double(x1) ;
        iter = iter + 1 ;
    end
    disp('Root is: ') ;
    disp(x1) ;
    figure(2)
    plot(err_vs_iter)
    title('Secant')
    xlabel('Iteration');
    ylabel('Error');
end
