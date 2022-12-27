function ff=testfunctions_2d(x,y,caseNum)

switch caseNum
    case 1 
        ff = @(x,y) 1-x+x; 
    case 2
        ff = @(x,y) x; 
    case 3
        ff = @(x,y) y; 
    case 4
        ff = @(x,y) x.^2;
    case 5
        ff = @(x,y) y.^2;
    case 6
        ff = @(x,y) x.*y;
    case 7 
        ff = @(x,y) (x+y)/2;
    case 8
        ff=@(x,y) x.^3+y.^4;
    case 9
        ff=@(x,y) x.^5+y.^6;
    case 10
        ff=@(x,y) (x-y).^9;
    case 11
        ff = @(x,y) x.^3./(1+y.^2);
    case 12
        ff = @(x,y) (1+3*x+4*y)./(1+x.*y);
    case 13
        ff =@(x,y) (1+x.^2+y.^2)./(5-x-y.^2);
    case 14
        ff =@(x,y) 1./(1+x.^2+y.^2);
    case 15
        ff = @(x,y)  (x.^10+y.^5)/2;
    case 16
        ff=@(x,y)  ((1+2*x+4*y)/7).^3;
    case 17
        ff= @(x,y) ((1+2*x+ 3*y)/6).^9;
    case 18 
        ff= @(x,y) ((1+2*x+ 3*y)/6).^5;
    case 19
        ff=@(x,y) sin(5*x+6*y);
    case 20
        ff=@(x,y) sin(2*pi*(x.^2+y.^2));
    case 21
        ff=@(x,y) tan(x-y/2);
    case 22
        ff=@(x,y) sin(x.^2-y.^2);
    case 23
        ff=@(x,y) cos(x.^2-y.^2);
    case 24
        ff=@(x,y) sin(1+cos(x.^2-y.^2));
    case 25
        ff=@(x,y) cos(1+sin(1+cos(x.^2-y.^2)));
    case 26
        ff=@(x,y) exp(-cos(1+sin(1+cos(x.^2-y.^2))));
    case 27
        ff=@(x,y) sin(sin(sin(sin(x.^2-y.^2))));
    case 28
        ff=@(x,y) sin(sin(sin(sin(sin(x.^2-y.^2)))));
    case 29
        ff=@(x,y) sin(sin(sin(sin(sin(sin(x.^2-y.^2))))));
    case 30
        ff=@(x,y) 1./(3+sin(sin(sin(sin(x.^2-y.^2)))));
    case 31
        ff=@(x,y) cos(5*x).*cos(4*y);
    case 32
        ff=@(x,y) sin(5*x).^2;
    case 33
        ff=@(x,y) tan(x).*tan(y);
    case 34
       ff=@(x,y) tan(x/2).*tan(y/2);
    case 35
        ff=@(x,y) (x-2*y).*sin(2*x-y);
    case 36
        ff=@(x,y) sin(10*x).*sin(10*y);
    case 37
        ff=@(x,y) (x.^2-y.^2)./(1+sin(x).^2+sin(y).^2);
    case 38
        ff=@(x,y) cos(x.^2-y.^2)./(1+x.^2+y.^3);
    case 39
        ff=@(x,y) sin(1+cos(x.^2-y.^2))./(1+cos(pi*x/2));
    case 40
        ff = @(x,y) sin(1+log(1+x.^2+y.^2)); 
    case 41 
        ff = @(x,y) log(1+x.^2+y.^2)/log(3); 
    case 42
        ff = @(x,y) x.*exp(1-x.^2-y.^2); 
    case 43
        ff = @(x,y) y.*log(1+x.^2+y.^2); 
    case 44
        ff = @(x,y) x.^2.*log(1+x.^2+y.^2);
    case 45
        ff = @(x,y) y.^2.*exp(1-x.^2-y.^2);
    case 46
        ff = @(x,y) x.*y.*log(1+x.^2+y.^2).*sin(pi*x+y);
    case 47 
        ff = @(x,y) sin(x+y).*log(1+x.^2+y.^2);
    case 48
        ff = @(x,y) log(x.^2-y.^3+2)/log(3);
    case 49
        ff = @(x,y) (log(5+ x.^3))./log(3+y.^2);
    case 50
        ff= @(x,y) log(1+x.^2+y.^2)./(1+x.^2+y.^2);
    case 51
        ff = @(x,y) (x+ y.^2).^2./(3+sin(x+y.^2)); 
    case 52
        ff = @(x,y) sin(x).^3./(3+sin(y).^2);
    case 53
        ff =@(x,y) sin(x+y)./(1+x.^2+y.^2);
    case 54
        ff = @(x,y) cos(x.^2+1./(1+x.*y));
    case 55
        ff=@(x,y) sin(1+2*x+3*y)./(2+sin(1+3*x+2*y));
    case 56
        ff =@(x,y) exp(1-x.^2-y.^2)./exp(1-x-y);
    case 57
        ff = @(x,y) (log(x.^2-y.^3+2)/log(3)).^2;
    case 58
        ff=@(x,y) cos(3*y).^2;
    case 59
        ff= @(x,y) ((1+2*x+ 3*y)/6).^7;
    case 60
        ff=@(x,y) sin(1+cos(x.^2+y.^2));
    case 61
        ff=@(x,y) (cos(1+sin(1+cos(x.^2+y.^2)))).^2;
    case 62
        ff=@(x,y) exp(-sin(1+sin(1+sin(x.^2+y.^2))));
    case 63
        ff=@(x,y) (sin(sin(sin(sin(sin(sin(x.^2+y.^2))))))).^3;
    case 64
        ff=@(x,y) sin(x+y)./(3+sin(sin(sin(sin(x.^2+y.^2)))));
    case 65    
        ff=@(x,y) sin(2*x.^2+3*y.^3);
    case 66
        ff = @(x,y) log(x.^2-y.^3+2)/log(3);
    case 67
        ff =@(x,y) exp(1-x.^2-y.^2)/exp(1);
    case 68
        ff=@(x,y) sin(x.^2+y.^2);
    case 69
        ff =@(x,y) sin(5*x).*sin(4*y);
    case 70
        ff=@(x,y) sin(2*x).*sin(3*y);
    case 71
        ff=@(x,y) sin(6*x).*sin(7*y);
    case 72
        ff=@(x,y) sin(8*x).*sin(8*y);
    case 73
        ff=@(x,y) cos(2*(x-y)/pi);
    case 74
        ff=@(x,y) max(x-0.5,0);
    case 75
        ff=@(x,y) max(x-0.5,0).*max(y-0.5,0);
    case 76
        ff=@(x,y) max(x-0.25,0) - max(x-0.75,0);
    case 77
        ff=@(x,y) (20*max(x-0.25,0).*max(0.75-x,0)).^2;
    case 78
        ff=@(x,y) (20*max(x-0.25,0).*max(0.75-x,0)).^2.*(20*max(y-0.25,0).*max(0.75-y,0)).^2;
    case 79
      ff = @(x,y) (100*(max(x-0.2,0).*max(0.8-x,0).*max(y-0.2,0).*max(0.8-y,0))).^2;
    case 80
      ff= @(x,y) (x./(1+exp(100*(x-0.25)))+ y./(1+exp(100*(0.75-x)))).* ...
          (1./(1+exp(100*(y-0.25)))+ 1./(1+exp(100*(0.75-y))));
    case 81
        ff=@(x,y) (max(x-0.25,0).*max(y-0.25,0)).^3;
    case 82
        ff = @(x,y) (max(x-0.25,0).*max(0.75-x,0)).^2;
    case 83
        ff=@(x,y) (max(x-0.3,0).*max(y-0.3,0)).^2;
    case 84
        ff=@(x,y)  tan(x.^2-y.^2)/tan(1);
    case 85
        ff=@(x,y)  (tan(x-y)).^2;
    case 86
        ff=@(x,y)  (tan(x-y)).^3;
    case 87
        ff=@(x,y)  exp(x.^3-y.^2);
    case 88
        ff=@(x,y)  sin(2*pi*(x+y));
    case 89
        ff=@(x,y) log(x.^2-y.^3+2);
    case 90
        ff=@(x,y) cos(1./(1+x.*y));
    case 91
        ff=@(x,y) cos(2*pi*(x.^2+y.^2));
    case 92
        ff=@(x,y) cos(2*pi*(x.^2-y.^2));
    case 93
        ff=@(x,y) sin(2*pi*x).*sin(2*pi*y);
    case 94
        ff=@(x,y) (x.*(1-x).*y.*(1-y)).^2;
    case 95
        ff=@(x,y) (sin(pi*x).*sin(pi*y)).^2;
    case 96
        ff=@(x,y) (x.*(1-x).*y.*(1-y));
    case 97 
        ff=@(x,y)  (x.*(1-x).*y.*(1-y)).^4;
    case 98
        ff=@(x,y) (sin(pi*x).*sin(pi*y)).^3;
    case 99
        ff=@(x,y) (sin(pi*x).*sin(pi*y)).^4;
    case 100
        ff=@(x,y)  (x.*(1-x).*y.*(1-y)).^8;
end