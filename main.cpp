//Code to solve the flow field vortex and compare the results to analytical solution:
//Written by: Prerna Patil
//Last edited: 12th May 2016
#define _USE_MATH_DEFINES
#include<cmath>
#include<stdio.h>
#include<iostream>
#include<math.h>
#include<iomanip>
#include "M4.h"
#include<fstream>
#include<iomanip>
#include<string.h>
using namespace std;
/**********************CLASS PARTICLE*************************************/
//Creating class Particle whcih stores information about X, Y, vorticity, u, v :
class particle
{
    double X, Y, vorticity, U, V, Gamma; //Gamma= Circulation
    public:
    //Default constructor:
    particle()
    {
        X=0.0; Y=0.0; vorticity=0.0; U=0.0; V=0.0; Gamma=0.0;
    }
    void set_position(double, double);
    void set_velocity(double, double);
    void set_vorticity(double);
    void set_Gamma(double);
    double getX();
    double getY();
    double getU();
    double getV();
    double getvor();
    double getGamma();
    //Operator overloading:
    void operator=(const particle& par)
    {
        this->X = par.X;
        this->Y = par.Y;
        this->U = par.U;
        this->V = par.V;
        this->vorticity = par.vorticity;
        this->Gamma = par.Gamma;
    }
};

double particle::getX(void)
{
    return X;
}

double particle::getY(void)
{
    return Y;
}

double particle::getU(void)
{
    return U;
}

double particle::getV(void)
{
    return V;
}

double particle::getvor(void)
{
    return vorticity;
}

double particle::getGamma(void)
{
    return Gamma;
}

void particle::set_position(double xin, double yin)
{
    X= xin;
    Y= yin;
}

void particle::set_velocity(double uin, double vin)
{
    U= uin;
    V= vin;
}

void particle::set_vorticity(double vor)
{
    vorticity =vor;
}

void particle::set_Gamma(double gam)
{
    Gamma= gam;
}
/********************** CLASS PARTICLE END *************************************/
int main()
{
    //Initialisation of the domain:
    //Initialisaing a Gaussian vortex field
    //\Sigma = 1 (centred at the origin)
    //Infinite domain??
    // x E [-3, 3] y E [-3, 3]
    double xstart = -5.0, xend = 5.0;
    double temp;
    double dt= 0.001, t=0.0, nu =0.1, epsilon = nu*dt; //nu= viscosity
    double T = 2.0; //Final time
    /*************************** INITIALIZATION **************************************/
    //Eulerian mesh
    int N=40;
    double h = (xend - xstart) / static_cast<double>(N-1);
    double x[N], y[N];
    //Volume remains unchanged for all time steps(incompressible flow)
    double Vol = pow((xstart - xend),2.0)/(N*N);
    for(int i=0; i<N; i++)
    {
        x[i] = xstart + i*h;
        y[i] = x[i] ;
        //cout<<"x: "<<x[i] <<endl;
    }
    clock_t begin, end;
    double time_spent;
    begin = clock();

    //Initialising particles on heap:
    //For using Adam-Bashforth second order for time stepping
    particle *par = new particle [N*N]; //Particle data at n+1 time step (current time step)
    particle *par_n = new particle [N*N]; //Particle data at n time step
    particle *par_n1 = new particle [N*N]; //Particle data at n-1 time step
    particle *par_diffuse = new particle[N*N]; //particles field data after diffusion
    /*
    //Initialise the reset mesh (all quatiies zero and set position)
    //This initialisation needs to be done before every reset
    for(int i=0; i<N; i++)
    {
    for(int j=0; j<N; j++)
    {
    par_reset[N*i+j].set_velocity(0.0, 0.0);
    par_reset[N*i+j].set_Gamma(0.0);
    par_reset[N*i+j].set_vorticity(0.0);
    par_reset[N*i+j].set_position(x[i], y[j]);
    }
    }
    */
    //Initialise the field at t=0;
    for(int i=0; i<N; i++)
    {
        for(int j=0; j<N; j++)
        {
            par_n1[N*i+j].set_position(x[i], y[j]);
            temp= 100.0* exp(-0.5*(pow(par_n1[N*i+j].getX(),2.0)+ pow(par_n1[N*i+j].getY(),2.0)));
            par_n1[N*i+j].set_vorticity(temp); //Setting the vorticity value for all the particles
            temp= par_n1[N*i+j].getvor() * Vol;
            par_n1[N*i+j].set_Gamma(temp); //Setting the Circulation value for all the particles
            //cout<<setw(2)<<par_n1[N*i+j].getGamma()<<" ";
        }
        //cout<<endl;
    }
    double kernal, r2, conX, conY, Circ;
    //We need to obtain the value of velocity at the grid points using the greens function:
    for(int p=0; p<N*N; p++)
    {
        conX = 0.0, conY = 0.0;
        for(int q=0; q< N*N; q++)
        {
            if(p==q){q=q+1;}
            r2 = pow(par_n1[q].getX() - par_n1[p].getX() ,2.0) + pow(par_n1[q].getY() - par_n1[p].getY() ,2.0);
            conX = conX - par_n1[q].getGamma() * (par_n1[p].getY() - par_n1[q].getY())/(2*M_PI*r2);
            conY = conY + par_n1[q].getGamma() * (par_n1[p].getX() - par_n1[q].getX())/(2*M_PI*r2);
        }
        par_n1[p].set_velocity(conX, conY);
    }
    //wrting initial data to file:
    ofstream outputtec;
    outputtec.open("tec.dat");
    outputtec<<"TITLE=\"Simple\""<<endl;
    outputtec<<"VARIABLES=\"X\" \"Y\" \"U\" \"V\" \"Gamma\""<<endl;
    outputtec<<"ZONE I="<<N<<" J="<<N<<endl;
    for(int i=0; i<N; i++)
    {
        for(int j=0; j<N; j++)
        {
            outputtec<<setw(6)<<par_n1[N*i+j].getX()<<" "<<setw(6)<<par_n1[N*i+j].getY()<<" "
            <<setw(6)<<par_n1[N*i+j].getU()<<" "<<setw(6)<<par_n1[N*i+j].getV()<<" "
            <<setw(6)<<par_n1[N*i+j].getGamma()<<endl;
        }
    }
    //Calculating the Gamma using the difussion equation: (for the first time step)
    /*
    //Convection by Adam bashforth second order:
    conX = par_n[p].getX() + dt* (1.5*par_n[p].getU() - 0.5*par_n1.getU());
    conY = par_n[p].getY() + dt* (1.5*par_n[p].getV() - 0.5*par_n1.getV());
    par[p].set_position(ConX, ConY);
    //Finding the velocities by MPE:
    r2 = pow(par_n1[q].getX() - par_n1[p].getX() ,2) + pow(par_n1[q].getY() - par_n1[p].getY() ,2);
    conX = -par[p].getGamma()*(par[q].getY() - par[p].getY())/(2*M_PI*r2);//Sum over p
    conY = par[p].getGamma() * (par[q].getX() - par[p].getX())/(2*M_PI*r2);//Sum over p
    par[p].set_velocity(conX, conY);
    */
    for(int p=0; p<N*N; p++)
    {
        Circ = par_n1[p].getGamma();
        conX = 0.0, conY = 0.0;
        for(int q=0; q<N*N; q++)
        {
            if(p==q){q=q+1;}
            r2 = pow(par_n1[q].getX() - par_n1[p].getX() ,2.0) + pow(par_n1[q].getY() - par_n1[p].getY() ,2.0);
            kernal = 1.0/(4.0*epsilon*M_PI) * exp(-r2/(4.0*epsilon));
            Circ = Circ + nu*dt * kernal * Vol*(par_n1[q].getGamma() - par_n1[p].getGamma());
            conX = conX - par_n1[q].getGamma() * (par_n1[p].getY() - par_n1[q].getY())/(2.0*M_PI*r2);
            conY = conY + par_n1[q].getGamma() * (par_n1[p].getX() - par_n1[q].getX())/(2.0*M_PI*r2);
        }

        par_diffuse[p].set_Gamma(Circ);
        par_diffuse[p].set_velocity(conX, conY);
        //Convection of the particles: (For the first time step: Euler forward)
        conX = par_n1[p].getX() + dt*par_n1[p].getU();
        conY = par_n1[p].getY() + dt*par_n1[p].getV();
        par_diffuse[p].set_position(conX, conY);
    }
    //Remeshing:
    //Remesh the diffused particles to the stable grid:
    //Interpolate the values of the circulation over the nearest grid points and regularize the grid:
    for(int i=0; i<N; i++)
    {
        for(int j=0; j<N; j++)
        {
            par_n[N*i+j].set_position(x[i], y[j]); //Setting the value of x and y of the new remesh grid
            //Since the first time step all values are initialised to zero by the default constructor
            //par_n[N*i+j].set_velocity(0.0, 0.0);
            //par_n[N*i+j].set_Gamma(0.0);
            //par_n[N*i+j].set_vorticity(0.0);
            temp =0.0;
            for(int p=0; p<N*N; p++)
            {
                //Using the M4 Third order accuracy
                temp = temp + par_diffuse[p].getGamma()* M4(par_n[N*i+j].getX() , par_diffuse[p].getX(), h) *
                M4(par_n[N*i+j].getY() , par_diffuse[p].getY(), h);
            }
            par_n[N*i+j].set_Gamma(temp); //Setting the value of the circulation at the 1st time step
            //cout<<setw(2)<<par_n[N*i+j].getGamma()<<" ";
        }
        //cout<<endl;
    }
    //We need to obtain the value of velocity at the grid points using the greens function:
    for(int p=0; p<N*N; p++)
    {
        conX = 0.0, conY = 0.0;
        for(int q=0; q< N*N; q++)
        {
            if(p==q){q=q+1;}
            r2 = pow(par_n[q].getX() - par_n[p].getX() ,2.0) + pow(par_n[q].getY() - par_n[p].getY() ,2.0);
            kernal = 1.0/(4.0*epsilon*M_PI) * exp(-r2/(4.0*epsilon));
            conX = conX - par_n[q].getGamma() * (par_n[p].getY() - par_n[q].getY())/(2.0*M_PI*r2);
            conY = conY + par_n[q].getGamma() * (par_n[p].getX() - par_n[q].getX())/(2.0*M_PI*r2);
        }
        par_n[p].set_velocity(conX, conY);
    }
    //Starting the time stepping using Adam Bashforth Second Order:
    t= dt;
    while(t < T)
    {
        //Initialise diffuse to zero:
        //Initialise the reset mesh (all quatiies zero and set position)

        for(int i=0; i<N; i++)
        {
            for(int j=0; j<N; j++)
            {
                par_diffuse[N*i+j].set_velocity(0.0, 0.0);
                par_diffuse[N*i+j].set_Gamma(0.0);
                par_diffuse[N*i+j].set_vorticity(0.0);
                par_diffuse[N*i+j].set_position(x[i], y[j]);
            }
        }
        //Diffuse the particles according to adam bashforth second order
        for(int p=0; p<N*N; p++)
        {
            Circ = par_n[p].getGamma();
            conX = 0.0, conY = 0.0;
            for(int q=0; q<N*N; q++)
            {
                if(p==q){q=q+1;}
                r2 = pow(par_n[q].getX() - par_n[p].getX() ,2.0) + pow(par_n[q].getY() - par_n[p].getY() ,2.0);
                kernal = 1.0/(4.0*epsilon*M_PI) * exp(-r2/(4.0*epsilon));
                Circ = Circ + nu*dt * kernal * Vol*(par_n[q].getGamma() - par_n[p].getGamma());
                conX = conX - par_n[q].getGamma() * (par_n[p].getY() - par_n[q].getY())/(2.0*M_PI*r2);
                conY = conY + par_n[q].getGamma() * (par_n[p].getX() - par_n[q].getX())/(2.0*M_PI*r2);
            }
            par_diffuse[p].set_Gamma(Circ);
            par_diffuse[p].set_velocity(conX, conY);
            //Convection of the particles: (Using Adam Bashforth 2nd order)
            conX = par_n[p].getX() + dt*(1.5*par_n[p].getU() - 0.5*par_n1[p].getU());
            conY = par_n[p].getY() + dt*(1.5*par_n[p].getV() - 0.5*par_n1[p].getV());
            par_diffuse[p].set_position(conX, conY);
            //cout<<par_diffuse[p].getU()<<" "<<endl;
        }
        //Remesh the particles to stable grid:
        for(int i=0; i<N; i++)
        {
            for(int j=0; j<N; j++)
            {
                par[N*i+j].set_position(x[i], y[j]); //Setting the value of x and y of the new remesh grid
                //Since the first time step all values are initialised to zero by the default constructor
                //par[N*i+j].set_velocity(0.0, 0.0);
                //par[N*i+j].set_Gamma(0.0);
                //par[N*i+j].set_vorticity(0.0);
                temp =0.0;
                for(int p=0; p<N*N; p++)
                {
                    //Using the M4 Third order accuracy
                    temp = temp + par_diffuse[p].getGamma()* M4(par[N*i+j].getX() , par_diffuse[p].getX(), h) *
                    M4(par[N*i+j].getY() , par_diffuse[p].getY(), h);
                }
                par[N*i+j].set_Gamma(temp); //Setting the value of the circulation at the 1st time step
                //cout<<par[N*i+j].getGamma()<<" ";
            }
            //cout<<endl;
        }
        //cout<<"***********"<<endl;
        //We need to obtain the value of velocity at the grid points using the greens function:
        for(int p=0; p<N*N; p++)
        {
            conX = 0.0, conY = 0.0;
            for(int q=0; q< N*N; q++)
            {
                if(p==q){q=q+1;}
                r2 = pow(par[q].getX() - par[p].getX() ,2.0) + pow(par[q].getY() - par[p].getY() ,2.0);
                kernal = 1.0/(4.0*epsilon*M_PI) * exp(-r2/(4.0*epsilon));
                conX = conX - par[q].getGamma() * (par[p].getY() - par[q].getY())/(2.0*M_PI*r2);
                conY = conY + par[q].getGamma() * (par[p].getX() - par[q].getX())/(2.0*M_PI*r2);
            }
            par[p].set_velocity(conX, conY);
        }
        //Shift the data for the particles:
        for(int p=0; p<N*N; p++)
        {
            par_n1[p] = par_n[p];
            par_n[p] = par[p];
        }
        t=t+dt; //Advance in the time step
        cout<<"t:"<<t<<endl;
        if(static_cast<int>(t/dt) %10 == 0)
        {
            //Write data to a file:
            //ofstream outputvor, outputU, outputV, outputGamma;
            //outputvor.open("Vorticity_"+ std::to_string(static_cast<int>(t/dt))+".dat");
            //outputU.open("U_"+ std::to_string(static_cast<int>(t/dt))+".dat");
            //outputV.open("V_"+ std::to_string(static_cast<int>(t/dt))+".dat");
            //outputGamma.open("Gamma_"+ std::to_string(static_cast<int>(t/dt))+".dat");
            outputtec<<"ZONE I="<<N<<" J="<<N<<endl;
            for(int i=0; i<N; i++)
            {
                for(int j=0; j<N; j++)
                {
                    // outputvor<<setw(6)<<par[N*i+j].getvor()<<" ";
                    // outputU<<setw(6)<<par[N*i+j].getU()<<" ";
                    // outputV<<setw(6)<<par[N*i+j].getV()<<" ";
                    outputtec<<setw(6)<<par[N*i+j].getX()<<" "<<setw(6)<<par[N*i+j].getY()<<" "
                    <<setw(6)<<par[N*i+j].getU()<<" "<<setw(6)<<par[N*i+j].getV()<<" "<<setw(6)<<par[N*i+j].getGamma()<<endl;
                    // outputGamma<<setw(6)<<par[N*i+j].getGamma()<<" ";
                }
                // outputvor<<endl;
                // outputU<<endl;
                // outputV<<endl;

                // outputGamma<<endl;
            }
        }
    }
    end = clock();
    time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    cout<<"Timetaken:"<<time_spent<<endl;
    return 0;
}