#include<iostream>
#include<vector>
//give ODE for solving
#define f(t,c)  c*(1-c) - s*c

using namespace std;

int main()

{
    //floats for initial conditions, parameters, etc.
    float t0, tfinal, c0, cn, dt, k1, k2, k3, k4, k, sgrid, sprobs;

    //ints for iterations/number of steps:
    int i, n;

    t0=0; //initial time
    tfinal=5; //final time
    c0=1; //initial condition
    n=5; //number of time steps

    //compute time step size dt:
    dt=(tfinal - t0)/n;

    //define normal dist
    double sigma, mu;
    sigma=.05;
    mu=.5;
    //(1/(sigma1*sqrt(2*pi)))*exp(-.5*((x - mu1)/sigma1).^2)
    double get_normal(sigma, mu, sgrid)
    {
        double normdist;
        normdist=(1/(sigma*sqrt(2*pi)))*exp(-.5*((sgrid - mu)/sigma).^2);
        return normdist;
    }

    //define s stuff
    vector<vector<double>>all_time;
    vector<vector<double>>all_space;

    //here comes big for loop
    for(int s = 0; s<12;++s)
    {
        double sgrid, sprobs;
        sgrid=s/12;
        sprobs=get_normal(mu,sigma,sgrid);

    //initialize time and space vecs
    vector<double>run_time;
    vector<double>run_space;

    //push back=append
    run_time.push_back(t0);
    run_space.push_back(c0);

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
        run_time.push_back(t0);
        run_space.push_back(c0);
    }
    all_time.push_back(run_time);
    all_space.push_back(run_space);
//big for loop ends here

    }
    //cout<<"\nValue of c at t=" <<tfinal <<" is " <<cn;

    //check size of run_time and run_space
    std::cout << "myvector stores " << int(all_space[3].size()) << " numbers.\n";

    return 0;

    //currently, this code takes in an inital and final t, as well as an initial c value
    //it spits out the final c value at the final t value

    //next step: apply true ode from research problem

    //next step: keep track of the whole solution curve at each time step
    //so need to figure out how to initalize vector and stuff

    //after that, need to loop over vector of sprobs and generate matrix of values

    //after that, need to choose a few different distributions to apply

    //after that, need to figure out how to time code in clion and in matlab

    //notes to self 11/22
    //need to fully define sgrid
    //need to figure out figures
    //need to get normal dist over sgrid
}