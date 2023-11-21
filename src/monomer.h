#ifndef MONOMER_H
#define MONOMER_H

#include <vector>

using namespace std;

class monomer
{
public:
    float x_pos, y_pos, z_pos;
    float x_force, y_force, z_force;
    vector<double> neighbors;

    monomer()
        : x_pos(0.0), y_pos(0.0), z_pos(0.0),
          x_force(0.0), y_force(0.0), z_force(0.0)
    {
    }

    vector<double> get_pos(), get_force();
    void reset_pos(), reset_force(), add_neighbor(double), print_neighbors();
};

void monomer::reset_pos()
{
    x_pos = 0.0;
    y_pos = 0.0;
    z_pos = 0.0;
}

void monomer::reset_force()
{
    x_force = 0.0;
    y_force = 0.0;
    z_force = 0.0;
}

vector<double> monomer::get_pos(){
    vector<double> result(3);
    result[0] = x_pos;
    result[1] = y_pos;
    result[2] = z_pos;
    return result;
}

vector<double> monomer::get_force(){
    vector<double> result(3);
    result[0] = x_force;
    result[1] = y_force;
    result[2] = z_force;
    return result;
}

void monomer::add_neighbor(double i){
    neighbors.push_back(i);
}

void monomer::print_neighbors(){
    for (int i = 0; i < neighbors.size(); i++){
        cout << neighbors[i] << " ";
    }
}
#endif