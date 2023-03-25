#include <algorithm>
#include <math.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <regex>
#include <vector>

/*--------------------------- Defining all variables ---------------------------*/

// Modifiables
const int event_count = 1000; // Sets maximum allowed choice for event_count

// Near modifiables
const int wwe_range = 1000; // Energy range of WW
const int wwp_range = 1000; // Z-Momentum range of W - 
float w1phi_com_selected = 0.0;

// Constants
const float invariant_mass = 13000; const float w_mass = 80.318; const float p_mass = 0.938272; // (M -> mass GeV/c^2)
// invariant mass ^ (GeV)

// Scalars
double px, py, pz, energy, mass; 
int pid, c1,c2,c3,c4,c5,c6,c7;
int direction_parameter; 
float wwpz_lab, wwe_lab;
float l1w1_angle_lab_candidate, l2w1_angle_lab_candidate, l1w2_angle_lab_candidate, l2w2_angle_lab_candidate;
float daw1, daw2;
float wp_candidate[6];
float wp_passing_candidate[6];
float w1px_com_selected, w1py_com_selected;
float w1p_selected;
float w1phi_lab, w2phi_lab;
float w1pxy_lab, w2pxy_lab;
float successRate;
int number_of_candidates;
float phi_increment;
int max_events;

double particle_lab[6][4][event_count];
double w1px_lab_event; double w1py_lab_event; double w1pz_lab_event; double w1e_lab_event;
double w2px_lab_event; double w2py_lab_event; double w2pz_lab_event; double w2e_lab_event;
double p1px_lab_event; double p1py_lab_event; double p1pz_lab_event; double p1e_lab_event;
double p2px_lab_event; double p2py_lab_event; double p2pz_lab_event; double p2e_lab_event;
double l1px_lab_event; double l1py_lab_event; double l1pz_lab_event; double l1e_lab_event;
double l2px_lab_event; double l2py_lab_event; double l2pz_lab_event; double l2e_lab_event;

float wwp_com;
float wwe_com;
float w1_gamma, w1_beta;
float w2_gamma, w2_beta;
float w1px_lab, w1py_lab, w1pz_lab; 
float w2px_lab, w2pz_lab, w2py_lab;
float w1e_lab, w2e_lab;
float wwe_diff;
float l1w1_angle_lab, l2w2_angle_lab;
string print_location;

// Vectors
std::vector<float> vector_w1phi_lab_candidate, vector_w2phi_lab_candidate;
std::vector<float> vector_w1w1, vector_w2w2, vector_w2w1, vector_w1w2;
std::vector<float> daw1_all,daw2_all;
std::vector<float> daw1_all_out,daw2_all_out,successRate_out,repeats,time_out;
std::vector<std::vector<float>> vector_wp_lab_candidate(6, std::vector<float>());
std::vector<std::vector<float>> wp_ref(3, std::vector<float>());

std::vector<float> sd1;
std::vector<float> sd2;
std::vector<float> mean1;
std::vector<float> mean2;

// ROOT Arrays
ROOT::Math::PxPyPzEVector w1_lab, w2_lab, p1_lab, p2_lab, l1_lab, l2_lab; 
ROOT::Math::PxPyPzEVector w1_lab_average, w2_lab_average;
ROOT::Math::PxPyPzEVector w1_com, w2_com, ww_com, w1_lab_ref, w2_lab_ref;
ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<double>,ROOT::Math::DefaultCoordinateSystemTag> w1_angle_diff, w2_angle_diff;
ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<double>,ROOT::Math::DefaultCoordinateSystemTag> w1_unit_distance_lab, w2_unit_distance_lab;
ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<double>,ROOT::Math::DefaultCoordinateSystemTag> w1_unit_distance_lab_ref, w2_unit_distance_lab_ref;
ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<double>,ROOT::Math::DefaultCoordinateSystemTag> l1_unit_distance_lab, l2_unit_distance_lab;
ROOT::Math::PxPyPzEVector w1_lab_candidate, w2_lab_candidate;
ROOT::Math::BoostZ wz_boost;
