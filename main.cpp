#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>
#include <complex>
#include "Matrix.h"
#include "parameters.h"
#include "rngesus.hpp"
#include "monomer.h"
#include <iomanip>
#include <cassert>
#include <string>
#include <sstream>
#include <cassert>

#define PI acos(-1.0)

using namespace std;

parameters prm;
xorshift64 rng;

double dot(vector<double> a, vector<double> b)
{
    double result;
    result = 0.0;
    for (int i = 0; i < 3; i++)
    {
        result += a[i] * b[i];
    }
    return result;
}

double norm(vector<double> a)
{
    return sqrt(dot(a, a));
}

vector<double> vecsum(vector<double> a, vector<double> b)
{
    vector<double> result(3);
    result[0] = a[0] + b[0];
    result[1] = a[1] + b[1];
    result[2] = a[2] + b[2];

    return result;
}

void printvec(vector<double> a)
{
    for (int i = 0; i < 3; i++)
    {
        cout << a[i] << " ";
    }
}

vector<double> vecneg(vector<double> a)
{
    vector<double> result(3);
    result[0] = -a[0];
    result[1] = -a[1];
    result[2] = -a[2];

    return result;
}

vector<double> displacement(monomer initial, monomer final) // initial to final vector
{
    return vecsum(final.get_pos(), vecneg(initial.get_pos()));
}

double distance(monomer a, monomer b)
{
    return norm(displacement(b, a));
}

vector<double> normalizevec(vector<double> v)
{
    vector<double> result(3);
    double n;
    n = norm(v);
    result[0] = (1.0 * result[0]) / n;
    result[1] = (1.0 * result[1]) / n;
    result[2] = (1.0 * result[2]) / n;
    return result;
}

vector<double> unitrandvec()
{
    vector<double> result(3), rvec(3);
    double nf;
    rvec[1] = rng.random();
    rvec[2] = rng.random();
    rvec[3] = rng.random();
    result = normalizevec(rvec);
    return result;
}

int main(int argc, char *argv[])
{
    if (argc != 2)
    {
        std::cerr << "Enter name of input file";
    }
    string inputfile = argv[1];
    prm.load(inputfile);

    // set seed
    rng.set_seed(prm.seed);

    vector<monomer> monomers(prm.length);

    // initialization

    for (int i = 0; i < monomers.size(); i++)
    {
        // reset saved values to zero
        monomers[i].reset_pos();
        monomers[i].reset_force();

        // set the initial position of each monomer along x-axis separated by bond_r0
        monomers[i].x_pos = i * prm.bond_r0;

        // initialize neighbors list
        if (i == 0)
        {
            monomers[i].add_neighbor(i + 1);
        }
        else if (i == monomers.size() - 1)
        {
            monomers[i].add_neighbor(i - 1);
        }
        else
        {
            monomers[i].add_neighbor(i - 1);
            monomers[i].add_neighbor(i + 1);
        }
    }

    ofstream dump;
    dump.open("dump.dat");

    double dx, dy, dz;
    double c1, c2;

    vector<double> r_in, r_i, r_n;
    double dist_in;

    double prefactor;

    vector<double> W(3);
    double Wnf;

    double sigma_by_r;

    c1 = prm.D / (prm.kB * prm.T);
    c2 = sqrt(2.0 * prm.D);

    for (int t_iter = 0; t_iter < prm.iterations; t_iter++)
    {
        // cout << "t = " << t_iter << endl;
        // At every iteration, reset forces (no inertia)
        for (int i = 0; i < monomers.size(); i++)
        {
            // cout << "Reset force on m" << i << endl;
            monomers[i].reset_force();
        }

        // Force simulation happens here

        for (int i = 0; i < monomers.size(); i++)
        {
            // cout << "Currently at m" << i << endl;

            // particle-particle forces
            for (int n = 0; n < monomers[i].neighbors.size(); n++)
            {

                // cout << "For monomer " << i << ", neighbor number " << n;
                // cout << " is monomer " << monomers[i].neighbors[n] << endl;

                r_i = monomers[i].get_pos();
                r_n = monomers[monomers[i].neighbors[n]].get_pos();
                r_in = vecsum(r_i, vecneg(r_n));
                dist_in = norm(r_in);

                if (prm.bond_break_message == 1)
                {
                    if (dist_in > 10.0 * prm.bond_r0)
                    {
                        cout << "Bond broken: monomers " << i << " and "
                             << monomers[i].neighbors[n]
                             << ", t = " << t_iter
                             << ", distance = "
                             << dist_in
                             << endl;
                    }
                }

                // spring force from harmonic bonds
                prefactor = -prm.bond_k * (dist_in - prm.bond_r0) * (1.0 / (dist_in));

                monomers[i].x_force += prefactor * r_in[0];
                monomers[i].y_force += prefactor * r_in[1];
                monomers[i].z_force += prefactor * r_in[2];

                // force from LJ-cut potential
                sigma_by_r = prm.LJ_sigma / dist_in;
                if (dist_in < prm.LJ_rc)
                {
                    prefactor = (48.0 / dist_in) * prm.LJ_e * (pow(sigma_by_r, 12) - 0.5 * pow(sigma_by_r, 6)) * (1.0 / (dist_in));

                    monomers[i].x_force += prefactor * r_in[0];
                    monomers[i].y_force += prefactor * r_in[1];
                    monomers[i].z_force += prefactor * r_in[2];
                }
            }

            // generate the noise
            W[0] = 2.0 * rng.random() - 1.0;
            W[1] = 2.0 * rng.random() - 1.0;
            W[2] = 2.0 * rng.random() - 1.0;

            Wnf = norm(W);

            W[0] = W[0] / Wnf;
            W[1] = W[1] / Wnf;
            W[2] = W[1] / Wnf;

            // calculate total displacements
            dx = c1 * monomers[i].x_force * prm.dt + c2 * W[0];
            dy = c1 * monomers[i].y_force * prm.dt + c2 * W[1];
            dz = c1 * monomers[i].z_force * prm.dt + c2 * W[2];

            // cap displacements

            // update positions
            monomers[i].x_pos += dx;
            monomers[i].y_pos += dy;
            monomers[i].z_pos += dz;
        }

        // Data dumping
        if (t_iter % prm.steps_to_skip == 0)
        {
            // cout << "Dumping at t = " << t_iter << endl;
            dump << prm.length << endl;
            dump << "polymer" << endl;

            for (int i = 0; i < monomers.size(); i++)
            {
                dump
                    << "1"
                    << " "
                    << monomers[i].x_pos << " "
                    << monomers[i].y_pos << " "
                    << monomers[i].z_pos << endl;
            }
        }
    }

    dump.close();

    return 0;
}