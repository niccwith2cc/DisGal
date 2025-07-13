# General Debugging in ```acoustic_solver_integrate ```

## $ M_{\text{inv}} $

I am not sure whether or not ```M_inv``` is the correct form in this case. ```M_loc``` consists of 2 series of the same general format.
Perhaps splitting it into two separate functions one for ```M_inv V``` and ``` M_inv P``` will make it better as we can also conduct it better

## Setting the initial condition

the way we are doing it now

```
v0=@(x) cos(pi*x);
p0=@(x) sin(pi*x);
w0 = [v0(x); p0(x)];
w = zeros(2*kp1*n, NT+1);
w(:, 1) = w0;
```

Is wrong to say the least. ``` w0 = [v0(x), p0(x)];``` should be the case first of all. Second of all

```
 w0 = [analytical_v(x,0); analytical_p(x,0)]
```
using this code would set the ```analytical_p ``` to always 0 which is correct as described by the initial conditions

## Runtime Loop

In this case, I think it is important to set up two different evolutions of p and v.

## Possible implementations

2 acoustic rhs equations, one for v and the other for p. Which would set up two different rhs integrations.

