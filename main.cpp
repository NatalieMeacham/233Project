#include<iostream>
#include<vector>
#include<cmath>
#include<numeric>
#include<fstream>
//give ODE for solving
#define f(t,c,s)  rho*c*(1-c) - d*s*c //included s to change its value in for loop later

using namespace std;

//declare normal dist function
void myFunction();
vector<double> Normal(vector<double> smesh, float sigma, float mu, float pi, double numpoints){
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
    //also create a text file to hold the probabilities vector
    ofstream probsfile;
    probsfile.open ("ProbVecFile.txt");
    double sprobsum2 = 0;
    int j;
    for (j = 0; j < numpoints; ++j) {
        sprobs[j] = sprobs[j] / sprobssum;
        //line to add value to text file
        //probsfile << "Writing this to a file.\n";
        //w_hid << std::to_string(W_hidden[a]) + "\n";
        probsfile<<std::to_string(sprobs[j]) + "\n";
        //cout<<"sprobs val here is "<<sprobs[i]<<endl;
        sprobsum2 = sprobsum2 + sprobs[j];
    }
    probsfile.close();
    cout<<"sprobsum2 is "<<sprobsum2<<endl;
    cout<<"New sprobs[8] is "<<sprobs[8]<<endl;
    cout<<"Size of sprobs is "<<sprobs.size()<<endl;

    return sprobs;

}

//linspace function to get vector smesh with numpoints number of points between x1 and x2.
vector<double> mylinspace(double x1, double x2, double numpoints){
    vector<double> smesh;
    double dgrid;
    smesh.assign(numpoints,0);
    int j;
    //dgrid = (x2 - x1)/numpoints;
    dgrid = (x2 - x1)/(numpoints-1); //try dgrid defined by matlab linspace
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
    //dt=(tfinal - t0)/n;
    dt=(tfinal - t0)/(n-1); //try the linspace formula
    //which one is correct?

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
    all_time.empty();
    all_space.empty();

    //define smesh (same thing as xgrid, different name in case it goes wrong)
    double numpoints, x1, x2;
    x1=0;
    x2=1;
    numpoints = 10;

    vector<double> smesh = mylinspace(x1, x2, numpoints); //create smesh

    //initialize time and space vecs
//    vector<double>run_time;
//    vector<double>run_space;


    //here comes big for loop over the number of s values
    for(int s = 0; s<numpoints;++s) //here we are looping over number of s values
    {
//    //initialize time and space vecs
    vector<double>run_time;
    vector<double>run_space;

        t0 = 0;
        c0 = 1;

//    //push back=append, apply to t0 and c0
//    run_time.push_back(t0);
//    run_space.push_back(c0);


    //RK method
    for (i=0; i<n; i++) //here we are looping over time for set s value
        //is i<n-1 more correct? error w matrix either way
    {
//        run_time.empty();
//        run_space.empty();

        //need to be specific about using specific s value
        //k1=dt*(f(t0,c0));
        //k2=dt*(f((t0 + dt/2), (c0 + k1/2)));
        //k3=dt*(f((t0 + h/2),(c0 + k2/2)));
        //k4=dt*(f((t0 + h),(c0 + k3)));
        //k=(k1 + 2*k2 + 2*k3 + k4)/6;
        //cn=c0 + k;
        //t0=t0 + dt; //reinitialize start time
        //c0=cn; //reinitialize start c val
        //run_time.push_back(t0);
        //

        run_time.push_back(t0);
        run_space.push_back(c0);

        k1=dt*(f(t0,c0,smesh[s]));
        k2=dt*(f((t0 + dt/2), (c0 + k1/2),smesh[s]));
        k3=dt*(f((t0 + h/2),(c0 + k2/2),smesh[s]));
        k4=dt*(f((t0 + h),(c0 + k3),smesh[s]));
        k=(k1 + 2*k2 + 2*k3 + k4)/6;
        cn=c0 + k;
        std::cout << "t0 = " << t0 << "\t c0 = " << c0 << "\t smesh = " << smesh[s] << "\t cn = " << cn << std::endl;
        t0=t0 + dt; //reinitialize start time
        c0=cn; //reinitialize start c val
        //run_time.push_back(t0);
        //run_space.push_back(cn);
    }
    all_time.push_back(run_time);
    all_space.push_back(run_space);
//big for loop ends here

    }
    //save matrix all_space to csv matrix
    ofstream outSol("MatSol.csv");
//    for (auto & row : all_space) {
//        for (auto col : row)
//            outSol << col <<',';
//        outSol << '\n';
//    }
////    for(int i = 0; i < n; i++){
//        for(int j =0; j < numpoints; j++){
//            outSol << all_space[j][i] << ',';
//        }
//        outSol << '\n';
//    }

    for(int i = 0; i < numpoints; i++){
        for(int j =0; j < n; j++){
            outSol << all_space[i][j] << ',';
        }
        outSol << '\n';
    }

ofstream outSolTime("MatSolTime.csv");
    for (auto & row : all_time) {
        for (auto col : row)
            outSolTime << col <<',';
        outSolTime << '\n';
    }

    //cout<<"\nValue of c at t=" <<tfinal <<" is " <<cn;

    //check size of run_time and run_space

    //how many entries are in 3rd row of all_space, aka how many time steps in solution for 3rd s value?
    std::cout << "myvector for space stores " << int(all_space[3].size()) << " numbers.\n";
    //std::cout << "myvector for space stores " << vector<double>(all_space.size()) << " numbers.\n";
    //std::cout << "myvector for time stores " << int(all_time[3].size()) << " numbers.\n";




    //once sprobs are calculated, can compute the weighted sum
    //need to extract rows of the matrix and then apply weights
vector<double> norm;
norm=Normal(smesh, sigma, mu, pi, numpoints); //apply normal distribution to smesh to get sprobs

    //initialize vector to hold final weighted solution
    vector<double> sumsol;
    sumsol.assign(n,0.);

    ofstream SumSol("SumSolFile.csv");

    for(int i =0; i < n; i++) {

        for (int s = 0; s < numpoints; ++s) //here we are looping over number of s values
        {
            sumsol[i] = sumsol[i] + all_space[s][i] * norm[s];


        }

        SumSol<< sumsol[i]<<'\n';
    }



        //try to get sample file to save to text file
    //ofstream mynewfile;
    //mynewfile.open ("hereisexample.txt");
    //mynewfile << "Here is some sample text.\n";
    //mynewfile.close();
    return 0;

    //return 0;

    //next step: keep track of the whole solution curve at each time step
    //so need to figure out how to initalize vector and stuff

    //after that, need to loop over vector of sprobs and generate matrix of values

    //after that, need to choose a few different distributions to apply

    //after that, need to figure out how to time code in clion and in matlab

    //notes to self 11/22
    //need to fully define sgrid
    //need to figure out figures
}



