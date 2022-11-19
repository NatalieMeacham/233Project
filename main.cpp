#include<iostream>

//give ODE for solving
#define f(t,c)  c*(1-c)

using namespace std;

int main()

{
    //floats for initial conditions, parameters, etc.
    float t0, tfinal, c0, cn, dt, k1, k2, k3, k4, k;

    //ints for iterations/number of steps:
    int i, n;

    t0=0; //initial time
    tfinal=5; //final time
    c0=1; //initial condition
    n=5; //number of time steps

    //compute time step size dt:
    dt=(tfinal - t0)/n;

    //RK method
    for (i=0; i<n; i++)
    {
        k1=dt*(f(t0,c0));
        k2=dt*(f((t0 + dt/2), (c0 + k1/2)));
        k3=dt*(f((t0 + h/2),(c0 + k2/2)));
        k4=dt*(f((t0 + h),(c0 + k3)));
        k=(k1 + 2*k2 + 2*k3 + k4)/6;
        cn=c0 + k;
        t0=t0 + dt; //reinitialize start time
        c0=cn; //reinitialize start c val
    }
    cout<<"\nValue of c at t=" <<tfinal <<" is " <<cn;

    return 0;

    //currently, this code takes in an inital and final t, as well as an initial c value
    //it spits out the final c value at the final t value

    //next step: apply true ode from research problem

    //next step: keep track of the whole solution curve at each time step
    //so need to figure out how to initalize vector and stuff

    //after that, need to loop over vector of sprobs and generate matrix of values

    //after that, need to choose a few different distributions to apply

    //after that, need to figure out how to time code in clion and in matlab
}