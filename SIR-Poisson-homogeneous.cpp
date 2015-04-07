/* Simulates independent realizations of a homogeneous SIR process on
a temporal network given on the form: (t i j) with one triple per line
and t,i, and j separated by tabs. ("\t").

Should be compiled with g++ including the boost library and using the
option -O2, i.e., as
g++ SIR-Poisson-homogeneous.cpp -o SIR -O2 -I<path of boost>

With the program compiled as SIR, it is called from the shell as:
./SIR <data> dt beta mu T_simulation ensembleSize outputTimeResolution
where:
<data> - path of text file containing contact data (temporal network);
dt - time-resolution of recorded contact data (time-step length);
beta - probability per time-step of infection when in contact with an
    infected node;
mu - probability per time-step of recovery for an infected node;
T_simulation - length of simulation in number of time-steps;
ensembleSize - number of independent realizations of the SIR process;
outputTimeResolution - time-resolution of the average number of infected
    and recovered nodes that the program gives as output.

The program gives as output three text files containing the average
number of infected nodes as function of time, the average number of
recovered nodes as function of time, and a histogram of the number of
recovered nodes at t=T_simulation-1.*/
//======================================================================
// Libraries
//======================================================================
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/exponential_distribution.hpp>
#include <algorithm> //for sort and unique
#include <vector>
#include <iterator>
#include <ctime>

//======================================================================
// Typedef
//======================================================================
// Variable types:
typedef unsigned int COUNTER;
typedef unsigned int NODE;
typedef std::vector<NODE> NODES; //list of nodes
typedef std::vector<bool> BOOLS;
struct CONTACT {NODE i; NODE j;}; //contact (i,j)
typedef std::vector<CONTACT> CONTACTS; //list of contacts during a single time-frame
typedef std::vector<CONTACTS> CONTACTS_LIST; //list of contact lists
// Random number generators:
typedef boost::mt19937_64 ENG; //use Mersenne Twister 19937 as PRNG engine
typedef boost::uniform_int<> DIST_INT; //define uniform distribution of integers
typedef boost::uniform_real<> DIST_REAL; //define uniform distribution of reals on [0,1)
typedef boost::exponential_distribution<> DIST_EXP; //define exponential distribution

//======================================================================
// Global parameters
//======================================================================
COUNTER N; //number of nodes in network
COUNTER T_data; //length of dataset (time steps)
COUNTER dt; //input time resolution (integer)
char inputname[200], outputname[200];
// Input/Output streams:
std::ofstream output;
std::ifstream input;
// Initialize PRNG:
ENG eng(9071982);
DIST_INT dist_int;
DIST_REAL dist_rand(0,1);
DIST_EXP dist_exp;

//======================================================================
// Function for importing list of contacts from 't\ti\tj\n' format file:
//======================================================================
CONTACTS_LIST loadContactListList(char *inputname)
{
    // Timer:
    std::clock_t clockStart = std::clock();

    // Define list of contact lists and list of nodes, etc.:
    std::string line;
    COUNTER t;
    NODE i,j;
    NODES nodes;
    CONTACTS_LIST contactListList;
    CONTACTS contactList;
    CONTACT contact;
    NODES::iterator last;

    // Print name of input file to screen:
    std::cout << "Filename: " << inputname << std::endl;

    input.open(inputname);
    // If input-file is not found or cannot be read, raise error and exit:
    if(input==NULL){ std::cout << "ERROR! File cannot be read."; return contactListList; }

    // Create list of nodes in file:
    while(getline(input,line))
    {
        input>>t>>i>>j;
        nodes.push_back(i);
        nodes.push_back(j);
    }
    T_data=t+1; //length of dataset
    std::cout << "T=" << T_data << std::endl;
    input.close();
    // Sort list and remove duplicates:
    std::sort(nodes.begin(),nodes.end());
    last=unique(nodes.begin(),nodes.end());
    nodes.resize(distance(nodes.begin(),last));
    N=nodes.size(); //number of unique nodes (network size)
    // List for redefining node IDs:
    NODES nodeIDs(nodes[nodes.size()-1]+1);
    std::cout << std::endl;
    for(COUNTER n=0; n<N; n++)
    {
        nodeIDs[nodes[n]]=n;
    }

    // Read first line of inputfile and get t,i,j:
    input.open(inputname);
    getline(input,line);
    input>>t>>i>>j;
    // Loop over t<T and create list of contact lists:
    for(COUNTER tt=0; tt<T_data; tt+=dt)
    {
        while(t==tt && !input.eof())
        {
            contact.i=nodeIDs[i]; contact.j=nodeIDs[j];
            contactList.push_back(contact);
            // Read line and get t,i,j:
            getline(input,line);
            input>>t>>i>>j;
        }
        contactListList.push_back(contactList);
        contactList.clear();
    }
    input.close();

    std::cout << std::endl << N << " nodes. Construction time: " << ( clock() - clockStart ) / (double) CLOCKS_PER_SEC << " s\n\n";

    return contactListList;
}

//======================================================================
// Main:
//======================================================================
int main(int argc, char *argv[])
{
    //-------------------------------------------------------------------------------------
    // Load contact data set:
    //-------------------------------------------------------------------------------------
    // Check if correct number of parameters was passed to program:
    if(argc<3){ std::cout << "Error! Data file not defined.\n"; return 0; }
    else
    {
        if(argc<5){ std::cout << "Error! Epidemic parameters not defined.\n"; return 0; }
        else
        {
            if(argc<6){ std::cout << "Error! Simulation time not specified.\n"; return 0; }
            else
            {
                if(argc<7){ std::cout << "Error! Ensemble size not specified.\n"; return 0; }
                else
                {
                    if(argc<8){ std::cout << "Error! Output time-resolution not specified.\n"; return 0;}
                }
            }
        }
    }
    // Set parameter values as specified:
    char *datafile=argv[1]; //dataset file name
    dt=atoi(argv[2]); //dataset time resolution (integer)
    char *beta_str = argv[3]; // base infection rate
    double beta = atof(beta_str);
    char *mu_str = argv[4]; // base recovery rate
    double mu = atof(mu_str);
    COUNTER betaprec=strlen(beta_str)-2;
    COUNTER muprec=strlen(mu_str)-2;
    COUNTER T_simulation = atoi(argv[5]); //simulation time
    COUNTER ensembleSize = atoi(argv[6]); //ensemble size (number of realizations)
    COUNTER outputTimeResolution = atoi(argv[7]); //output time-resolution

    // Open input file and load contact_lists:
    sprintf(inputname,"%s",datafile);
    CONTACTS_LIST contactListList=loadContactListList(inputname);
    // Check if length of contacts_list > 0 and end program if it is not:
    if(contactListList.size()==0){ std::cout << "Error! Dataset empty.\n"; return 0; }

    //-------------------------------------------------------------------------------------
    // Define variables:
    //-------------------------------------------------------------------------------------
    NODES infected; //list of infected nodes
    double Mu; //cumulative recovery rate
    BOOLS isSusceptible, isInfected; //list which nodes are susceptible/infected, respectively
    COUNTER I,R; //number of infected and recovered nodes
    COUNTER SI; //number of susceptible nodes in contact with infectious nodes
    NODES si_s; //list of susceptible nodes in contact with infected nodes
    double Beta; //cumulative infection rate
    double Lambda; //cumulative transition rate
    double xi; //
    COUNTER t; //time counter
    COUNTER t_infectionStart; //starting time of infection
    NODE root; //root node of infection
    NODE i,j;
    double tau; //renormalized waiting time until next event
    CONTACTS_LIST::iterator contactList_iterator; //iterator over list of contacts
    CONTACTS::iterator contact_iterator; //iterator over contacts
    double r_transitionType; //random variable for choosing which type of transition happens
    COUNTER m; //transition process number
    COUNTER n; //time-counter
    NODES::iterator node_iterator; //iterator over list of nodes
    NODES::iterator last; //iterator for use when generating unique list of new infected nodes
    // Containers for output data:
    NODES sumI_t(T_simulation/outputTimeResolution); //list of number of infected nodes in each recorded frame
    NODES sumR_t(T_simulation/outputTimeResolution); //list of number of recovered nodes in each recorded frame
    NODES hist_R(N+1); //histogram of R values after I=0
    // Random number generators:
    boost::variate_generator<ENG,DIST_INT> randint(eng,dist_int); //random integer
    boost::variate_generator<ENG,DIST_REAL> rand(eng,dist_rand); //random float on [0,1[
    boost::variate_generator<ENG,DIST_EXP> randexp(eng,dist_exp); //random exponentially distributed float

    //-------------------------------------------------------------------------------------
    // Simulate:
    //-------------------------------------------------------------------------------------
    std::clock_t clockStart = std::clock();     //timer
    COUNTER stopped=0; //counter of number of simulations that stopped (I=0) during T_simulation
    for(COUNTER q=0; q<ensembleSize; q++)
    {
        std::cout << q << "/" << ensembleSize << std::endl; //print realization # to screen
        // Choose at random infectious root node:
        root=randint(N);
        // Initialize lists of infected nodes and infected node IDs:
        R=0;
        infected.clear();
        infected.push_back(root);
        I=1;
        Mu=mu;
        isInfected.assign(N,false);
        isInfected[root]=true;
        isSusceptible.assign(N,true);
        isSusceptible[root]=false;
        // First waiting "time":
        tau=randexp(1);
        // Random starting time of infection:
        t_infectionStart=randint(T_data);
        // Reset simulation time to zero:
        t=0;

        // Loop until either I=0 or t>=T_simulation
        while(I>0 && t<T_simulation)
        {
            // Loop over list of contact lists:
            for(contactList_iterator=contactListList.begin()+t_infectionStart; contactList_iterator!=contactListList.end(); contactList_iterator++)
            {
                // Create list of susceptible nodes in contact with infected nodes:
                si_s.clear();
                for(contact_iterator=(*contactList_iterator).begin(); contact_iterator!=(*contactList_iterator).end(); contact_iterator++)
                {
                    i=(*contact_iterator).i;
                    j=(*contact_iterator).j;
                    if(isInfected[i])
                    {
                        if(isSusceptible[j])
                        {
                            si_s.push_back(j);
                        }
                    }
                    else
                    {
                        if(isInfected[j])
                        {
                            if(isSusceptible[i])
                            {
                                si_s.push_back(i);
                            }
                        }
                    }
                }
                SI=si_s.size(); //number of possible S->I transitions
                Beta=(double)SI*beta; //cumulative infection rate
                Lambda=Beta+Mu; //cumulative transition rate:

                // Check if transition takes place during time-step:
                if(tau>=Lambda) //no transition takes place
                {
                    tau-=Lambda;
                }
                else //at least one transition takes place
                {
                    xi=1.; //fraction of time-step left before transition
                    // Sampling step:
                    while(tau<xi*Lambda) //repeat if next tau is smaller than ~ Lambda-tau
                    {
                        xi-=tau/Lambda; //fraction of time-step left after transition
                        r_transitionType=Lambda*rand(); //random variable for weighted sampling of transitions
                        if(r_transitionType<Beta) //S->I
                        {
                            m=randint(SI); //transition m
                            isInfected[si_s[m]]=true;
                            isSusceptible[si_s[m]]=false;
                            // Add infected node to list:
                            infected.push_back(si_s[m]);
                            I++;
                            Mu+=mu;
                        }
                        else //I->R
                        {
                            m=randint(I); //transition m
                            isInfected[infected[m]]=false;
                            // Remove drawn element from infected:
                            infected[m]=infected.back();
                            infected.pop_back();
                            I--;
                            R++;
                            Mu-=mu;
                        }
                        // Redo list of S-I contacts:
                        si_s.clear();
                        for(contact_iterator=(*contactList_iterator).begin(); contact_iterator!=(*contactList_iterator).end(); contact_iterator++)
                        {
                           i=(*contact_iterator).i;
                           j=(*contact_iterator).j;
                           if(isInfected[i])
                           {
                                if(isSusceptible[j])
                                {
                                    si_s.push_back(j);
                                }
                            }
                            else
                            {
                                if(isInfected[j])
                                {
                                    if(isSusceptible[i])
                                    {
                                        si_s.push_back(i);
                                    }
                                }
                            }
                        }
                        SI=si_s.size();
                        Beta=beta*(double)SI;
                        Lambda=Beta+Mu; //new cumulative transition rate
                        // Draw new renormalized waiting time:
                        tau=randexp(1);
                    }
                }
                // Stop if I=0:
                if(I==0)
                {
                    stopped++;
                    hist_R[R]++;
                    for(n=t; n<T_simulation; n++)
                    {
                        if(n % outputTimeResolution ==0){ sumR_t[n/outputTimeResolution]+=R; }
                    }
                   break;
                }
                // Read out I and R if t is divisible by outputTimeResolution
                if(t % outputTimeResolution ==0)
                {
                    if(t>=T_simulation) //stop if max simulation time-steps has been reached
                    {
                        break;
                    }
                    else
                    {
                        sumI_t[t/outputTimeResolution]+=I;
                        sumR_t[t/outputTimeResolution]+=R;
                    }
                }
                t++;
            }
            t_infectionStart=0;
        }
    }

    //-------------------------------------------------------------------------------------
    // Save epidemic data to disk:
    //-------------------------------------------------------------------------------------
    double t_simu = ( std::clock() - clockStart ) / (double) CLOCKS_PER_SEC;

    clockStart=std::clock();
    // Open output file:
    sprintf(outputname,"sum(I_t)-%s,N=%u,dt=%u,T=%u,beta=%.*f,mu=%.*f,Q=%u,res=%u.txt",datafile,N,dt,T_simulation,betaprec,beta,muprec,mu,ensembleSize,outputTimeResolution);
    output.open(outputname);
    // Write avg(I_t) to file:
    for(node_iterator=sumI_t.begin(); node_iterator!=sumI_t.end(); node_iterator++)
    {
        output << *node_iterator << "\t";
    }
    output.close();

    // Open output file:
    sprintf(outputname,"sum(R_t)-%s,N=%u,dt=%u,T=%u,beta=%.*f,mu=%.*f,Q=%u,res=%u.txt",datafile,N,dt,T_simulation,betaprec,beta,muprec,mu,ensembleSize,outputTimeResolution);
    output.open(outputname);
    // Write avg(R_t) to file:
    for(node_iterator=sumR_t.begin(); node_iterator!=sumR_t.end(); node_iterator++)
    {
        output << *node_iterator << "\t";
    }
    output << "\n";
    output.close();
    // Open output file:
    sprintf(outputname,"h(R)-%s,N=%u,dt=%u,T=%u,beta=%.*f,mu=%.*f,Q=%u,res=%u.txt",datafile,N,dt,T_simulation,betaprec,beta,muprec,mu,ensembleSize,outputTimeResolution);
    output.open(outputname);
    // Write histogram of R, p(R) to file:
    for(node_iterator=hist_R.begin(); node_iterator!=hist_R.end(); node_iterator++)
    {
        output << *node_iterator << "\t";
    }
    output.close();
    double t_write = ( std::clock() - clockStart ) / (double) CLOCKS_PER_SEC;

    std::cout << std::endl << "temporal Gillespie---homogeneous & Poissonian SIR: N=" << N << ", T=" << T_data << ", beta=" << beta << ", mu=" << mu;
    std::cout << ", output time-resolution = " << outputTimeResolution << std::endl;
    std::cout << "Simulation time: " << t_simu << "s, Stopped simulations: " << stopped << "/" << ensembleSize << std::endl;
    std::cout << "Writing to file: " << t_write << "s" << std::endl;

    return 0;
}
