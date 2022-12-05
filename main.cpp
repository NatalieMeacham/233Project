#include<iostream>
#include<vector>
#include<cmath>
#include<numeric>
//give ODE for solving
#define f(t,c)  rho*c*(1-c) - d*s*c

using namespace std;

//declare normal dist function
void myFunction();
void Normal(vector<double> smesh, float sigma, float mu, float pi, double numpoints){
    vector<double> sprobs; //need to set it up as a vector to be able to take the size
    sprobs.assign(numpoints,0.);
    int i;
    double sprobssum = 0;
    for (i = 0; i < numpoints; ++i) {
        //takes in smesh and spits out normal equation as sprobs
        sprobs[i] = (1/(sigma*sqrt(2*pi)))*exp(-.5*((smesh[i] - mu)/sigma)*((smesh[i] - mu)/sigma));
        //cout<<"sprobs val here is "<<sprobs[i]<<endl;
        sprobssum = sprobssum + sprobs[i];
    }
    cout<<"sprobssum is "<<sprobssum<<endl;
    cout<<"Old sprobs[8] is "<<sprobs[8]<<endl;

    //normalize probabilities with sprobs[i] = sprobs[i]/sum(sprobs)

    double sprobsum2 = 0;
    int j;
    for (j = 0; j < numpoints; ++j) {
        sprobs[j] = sprobs[j] / sprobssum;
        //cout<<"sprobs val here is "<<sprobs[i]<<endl;
        sprobsum2 = sprobsum2 + sprobs[j];
    }
    cout<<"sprobsum2 is "<<sprobsum2<<endl;
    cout<<"New sprobs[8] is "<<sprobs[8]<<endl;
    cout<<"Size of sprobs is "<<sprobs.size()<<endl;

}

//linspace function to get vector smesh with numpoints number of points between x1 and x2.
vector<double> mylinspace(double x1, double x2, double numpoints){
    vector<double> smesh;
    double dgrid;
    smesh.assign(numpoints,0);
    int j;
    dgrid = (x2 - x1)/numpoints;
    for (j=0; j<numpoints; ++j){
        smesh[j] = x1 + dgrid*j;
    }
    cout<< "Scalar dgrid is "<< dgrid<<endl;
    cout<< "Smesh 3rd grid point is " <<smesh[3]<<endl;

    return smesh;
}


int main()

{
    //floats for initial conditions, parameters, etc.
    float t0, tfinal, c0, cn, dt, k1, k2, k3, k4, k, rho, d, mu, sigma, pi;

    //ints for iterations/number of steps:
    int i, n;

    t0=0; //initial time
    tfinal=5; //final time
    c0=1; //initial condition
    n=5; //number of time steps

    //compute time step size dt:
    dt=(tfinal - t0)/n;

    //parameters for ODE:
    rho=1;
    d=0.9;

    //parameters for normal distribution:
    mu=0.5;
    sigma=0.05;
    pi=3.14159;

    //define s stuff: time and space vectors?
    vector<vector<double>>all_time;
    vector<vector<double>>all_space;

    //define smesh (same thing as xgrid, different name in case it goes wrong)
    double numpoints, x1, x2;
    x1=0;
    x2=1;
    numpoints = 10;

    vector<double> smesh = mylinspace(x1, x2, numpoints); //create smesh

    //here comes big for loop over the number of s values
    for(int s = 0; s<numpoints;++s) //here we are looping over number of s values
    {
    //initialize time and space vecs
    vector<double>run_time;
    vector<double>run_space;

    //push back=append, apply to t0 and c0
    run_time.push_back(t0);
    run_space.push_back(c0);

    //RK method
    for (i=0; i<n; i++) //here we are looping over time for set s value
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

    //how many entries are in 3rd row of all_space, aka how many time steps in solution for 3rd s value?
    std::cout << "myvector for space stores " << int(all_space[3].size()) << " numbers.\n";
    //std::cout << "myvector for time stores " << int(all_time[3].size()) << " numbers.\n";




    //once sprobs are calculated, can compute the weighted sum
    //need to extract rows of the matrix and then apply weights

    Normal(smesh, sigma, mu, pi, numpoints); //apply normal distribution to smesh to get sprobs

    return 0;

    //next step: keep track of the whole solution curve at each time step
    //so need to figure out how to initalize vector and stuff

    //after that, need to loop over vector of sprobs and generate matrix of values

    //after that, need to choose a few different distributions to apply

    //after that, need to figure out how to time code in clion and in matlab

    //notes to self 11/22
    //need to fully define sgrid
    //need to figure out figures
}

// Function definition
//void myFunction() {
    //this function needs to take in an array of s values, called sgrid, and apply normal dist.

   // cout << "I just got executed!";
//}

