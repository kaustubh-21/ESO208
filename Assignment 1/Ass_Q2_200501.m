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
    disp('1. Muller Method');
    disp('2. Bairstow Method');
    
method_select = input('Select the method for finding Roots: ') ;
%Muller = 1, Bairstow = 2

if(method_select == 1)   %Muller Method
    x0 = input('Enter the x0 :  ');
    x1 = input('Enter the x1 :  ');
    x2 = input('Enter the x2 :  ');
    max_err = input('Enter maximum Error: ') ;
    max_iter = input('Enter maximum Iteration allowed: ') ;
    error = 100 ;
    iter = 0 ;
    while ((error > max_err) && (iter < max_iter))
        diff1 = (func(x1)-func(x0))/(x1-x0) ;
        diff2 = (func(x2)-func(x1))/(x2-x1) ;
        a = (diff2-diff1)/(x2-x0) ;
        b = diff2 + (x2-x1)*a ;
        c = func(x2) ;
        
        x3 = x2 - 2*c/(b + sign(b)*sqrt(b^2 - 4*a*c)) ;
        
        error = abs(100*(x3-x2)/x3) ;
        err_vs_iter(iter+1) = error ;
        x0 = x1 ;
        x1 = x2 ;
        x2 = x3 ;
        iter = iter + 1 ;
    end
    
    figure(2)
    plot(err_vs_iter) 
    title('Muller')
    xlabel('Iteration');
    ylabel('Error');
    disp('Root is : ') ;
    disp(x2) ;
end

if(method_select == 2)   %Bairstow Method
    
end