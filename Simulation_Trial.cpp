#include <algorithm>
#include <math.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <regex>
#include <vector>
#include <chrono>
#include "Simulation_Head_Reduced.h"
#include "Simulation_Style_Reduced.c"

void print(std::vector <float> &a) {
    std::cout << "The vector elements are : ";

    for(int i=0; i < a.size(); i++)
    std::cout << a.at(i) << ' ';
    }

int load_data(){ 
    std::ifstream infile;
    infile.open("Simulation_Data.txt"); 
    std::string line;       
    for(int i = 0; i < event_count; i++){
        std::getline(infile, line); 
        for(int j = 0; j < 10; j++){
            std::getline(infile, line);
            std::istringstream iss(line);
            iss >> pid >> c1 >> c2 >> c3 >> c4 >> c5 >> px >> py >> pz >> energy >> mass >> c6 >> c7;
            if(j==1&&pid==2212){particle_lab[0][0][i] = px; particle_lab[0][1][i] = py; particle_lab[0][2][i] = pz; particle_lab[0][3][i] = energy;}
            if(j==2&&pid==2212){particle_lab[1][0][i] = px; particle_lab[1][1][i] = py; particle_lab[1][2][i] = pz; particle_lab[1][3][i] = energy;} 
            if(pid==24){particle_lab[2][0][i] = px; particle_lab[2][1][i] = py; particle_lab[2][2][i] = pz; particle_lab[2][3][i] = energy;} // Strictly for comparison
            if(pid==-24){particle_lab[3][0][i] = px; particle_lab[3][1][i] = py; particle_lab[3][2][i] = pz; particle_lab[3][3][i] = energy;} // Strictly for comparison
            if(pid==13){particle_lab[4][0][i] = px; particle_lab[4][1][i] = py; particle_lab[4][2][i] = pz; particle_lab[4][3][i] = energy;} 
            if(pid==-13){particle_lab[5][0][i] = px; particle_lab[5][1][i] = py; particle_lab[5][2][i] = pz; particle_lab[5][3][i] = energy;}
            }//j
        std::getline(infile, line);
        }//i
    infile.close();

    return 0;
    }//load_data()

int create_directory(TTree * t){

    t->Branch("daw1", &daw1); // Decay opening angle of W +
    t->Branch("daw2", &daw2); // Decay opening angle of W -
    t->Branch("w1px_candidate", &wp_candidate[0]); // X-Momentum of W +
    t->Branch("w1py_candidate", &wp_candidate[1]); // X-Momentum of W -
    t->Branch("w1pz_candidate", &wp_candidate[2]); // Y-Momentum of W +
    t->Branch("w2px_candidate", &wp_candidate[3]); // Y-Momentum of W -
    t->Branch("w2py_candidate", &wp_candidate[4]); // Z-Momentum of W +
    t->Branch("w2pz_candidate", &wp_candidate[5]); // Z-Momentum of W -

    return 0;
    }//create_directory()

int reference_variables(int event){ 
    
    p1px_lab_event = particle_lab[0][0][event]; p1py_lab_event = particle_lab[0][1][event]; p1pz_lab_event = particle_lab[0][2][event]; p1e_lab_event = particle_lab[0][3][event];
    p2px_lab_event = particle_lab[1][0][event]; p2py_lab_event = particle_lab[1][1][event]; p2pz_lab_event = particle_lab[1][2][event]; p2e_lab_event = particle_lab[1][3][event];
    w1px_lab_event = particle_lab[2][0][event]; w1py_lab_event = particle_lab[2][1][event]; w1pz_lab_event = particle_lab[2][2][event]; w1e_lab_event = particle_lab[2][3][event];
    w2px_lab_event = particle_lab[3][0][event]; w2py_lab_event = particle_lab[3][1][event]; w2pz_lab_event = particle_lab[3][2][event]; w2e_lab_event = particle_lab[3][3][event];
    l1px_lab_event = particle_lab[4][0][event]; l1py_lab_event = particle_lab[4][1][event]; l1pz_lab_event = particle_lab[4][2][event]; l1e_lab_event = particle_lab[4][3][event];
    l2px_lab_event = particle_lab[5][0][event]; l2py_lab_event = particle_lab[5][1][event]; l2pz_lab_event = particle_lab[5][2][event]; l2e_lab_event = particle_lab[5][3][event];

    wp_ref[0].push_back(w1px_lab_event);
    wp_ref[1].push_back(w1py_lab_event);
    wp_ref[2].push_back(w1pz_lab_event);

    return 0;
    }//reference_variables()

int fill_lab_values(){ // With references

    // Comparison variables
    w1_lab_ref.SetCoordinates(w1px_lab_event, w1py_lab_event, w1pz_lab_event, w1e_lab_event);
	w2_lab_ref.SetCoordinates(w2px_lab_event, w2py_lab_event, w2pz_lab_event, w2e_lab_event);
    w1_unit_distance_lab_ref = w1_lab_ref.Vect().Unit(); // Unit direction vector based of momentum cop_massonents -> Direction of W +
    w2_unit_distance_lab_ref = w2_lab_ref.Vect().Unit(); // Unit direction vector based of momentum cop_massonents -> Direction of W -

    // Reconstruction variables
    p1_lab.SetCoordinates(p1px_lab_event, p1py_lab_event, p1pz_lab_event, p1e_lab_event);
    p2_lab.SetCoordinates(p2px_lab_event, p2py_lab_event, p2pz_lab_event, p2e_lab_event);
	l1_lab.SetCoordinates(l1px_lab_event, l1py_lab_event, l1pz_lab_event, l1e_lab_event);
	l2_lab.SetCoordinates(l2px_lab_event, l2py_lab_event, l2pz_lab_event, l2e_lab_event);
	l1_unit_distance_lab = l1_lab.Vect().Unit(); // Unit direction vector based of momentum cop_massonents -> Direction of Lepton 1
    l2_unit_distance_lab = l2_lab.Vect().Unit(); // Unit direction vector based of momentum cop_massonents -> Direction of Lepton 2

    return 0;
    }//fill_lab_values()

int add_smearing(float smearing_extent){

    float sigma_p1e = 2.5 + ((0.5*invariant_mass-p1e_lab_event)/p1e_lab_event) * 50;
    float sigma_p2e = 2.5 + ((0.5*invariant_mass-p2e_lab_event)/p2e_lab_event) * 50;

    sigma_p1e = p1e_lab_event + smearing_extent * gRandom->Gaus(0.0,sigma_p1e);
    sigma_p2e = p2e_lab_event + smearing_extent * gRandom->Gaus(0.0,sigma_p2e);
    
    wwe_lab = invariant_mass - sigma_p1e - sigma_p2e;
    p1_lab.SetCoordinates(p1px_lab_event, p1py_lab_event, p1pz_lab_event, sigma_p1e);//smearing
    p2_lab.SetCoordinates(p2px_lab_event, p2py_lab_event, p2pz_lab_event, sigma_p2e);//smearing

    return 0;
    }//add_smearing(smearing_extent)

int ww_boost_to_com(float smearing_extent){

    // Reconstructed
    wwpz_lab = p1pz_lab_event + p2pz_lab_event; // Combined Z-Momentum of protons must be equal to combined WW Z-Momentum
    if(wwpz_lab<=0.0){direction_parameter = -1;}
    else{direction_parameter = 1;}
    wwe_lab = invariant_mass - (p1e_lab_event + p2e_lab_event); // Energy of WW 

    //add_smearing(smearing_extent); //UNCOMMENT IF WANT TO USE

    ww_com.SetCoordinates(-(p1px_lab_event+p2px_lab_event), -(p1py_lab_event+p2py_lab_event), -wwpz_lab, wwe_lab);
    wz_boost.SetBeta(-direction_parameter*ww_com.Beta()); // BoostZ( beta:VALUE , gamma:VALUE )

    return 0;
    }//ww_boost_to_com(smearing_extent)

int w1w2_boost_to_lab(){ 

    // To compare: From the given data
    w1p_selected = sqrt(wwe_com*wwe_com/4-w_mass*w_mass - wwp_com*wwp_com);
    w1px_com_selected = w1p_selected*cos(w1phi_com_selected); // w1phi_com_selected set to 0 in header
    w1py_com_selected = w1p_selected*sin(w1phi_com_selected); // w1phi_com_selected set to 0 in header
    w1_com.SetCoordinates( w1px_com_selected,  w1py_com_selected,  wwp_com, wwe_com/2);
    w2_com.SetCoordinates(-w1px_com_selected, -w1py_com_selected, -wwp_com, wwe_com/2); // From symmetry

    // Can you compress into a 3-dimensional matrix
    w1_lab = wz_boost(w1_com);
    w2_lab = wz_boost(w2_com);
	w1e_lab = w1_lab.E(); w1px_lab = w1_lab.Px(); w1py_lab = w1_lab.Py(); w1pz_lab = w1_lab.Pz();
    w2e_lab = w2_lab.E(); w2px_lab = w2_lab.Px(); w2py_lab = w2_lab.Py(); w2pz_lab = w2_lab.Pz();

    w1_beta = sqrt(w1px_lab*w1px_lab + w1py_lab*w1py_lab + w1pz_lab*w1pz_lab)/w1e_lab; 
    w2_beta = sqrt(w2px_lab*w2px_lab + w2py_lab*w2py_lab + w2pz_lab*w2pz_lab)/w2e_lab;
    w1_gamma = 1/sqrt(1-w1_beta*w1_beta);
    w2_gamma = 1/sqrt(1-w2_beta*w2_beta);
	l1w1_angle_lab = acos( 1/w1_beta - w1_gamma * w_mass*(1/w1_beta - w1_beta)/2/l1e_lab_event ); 
    l2w2_angle_lab = acos( 1/w2_beta - w2_gamma * w_mass*(1/w2_beta - w2_beta)/2/l2e_lab_event ); 
			
    return 0;
    }//w1w2_boost_to_com()

int get_angle_candidates(float phi_increment){

    // Reconstructed                
    w1phi_lab = w1_lab.Phi();
    w1pxy_lab = sqrt(w1px_lab*w1px_lab + w1py_lab*w1py_lab);
    w1_lab_candidate.SetCoordinates(w1pxy_lab*cos(w1phi_lab+phi_increment), w1pxy_lab*sin(w1phi_lab+phi_increment), w1pz_lab, w1e_lab); // Rotated four-momenta
    w1_unit_distance_lab = w1_lab_candidate.Vect().Unit();
    l1w1_angle_lab_candidate = acos(l1_unit_distance_lab.Dot(w1_unit_distance_lab));

    w2phi_lab = w2_lab.Phi();
    w2pxy_lab = sqrt(w2px_lab*w2px_lab + w2py_lab*w2py_lab);
    w2_lab_candidate.SetCoordinates(w2pxy_lab*cos(w2phi_lab+phi_increment), w2pxy_lab*sin(w2phi_lab+phi_increment), w2pz_lab, w2e_lab); // Rotated four-momenta
    w2_unit_distance_lab = w2_lab_candidate.Vect().Unit();
    l2w2_angle_lab_candidate = acos(l2_unit_distance_lab.Dot(w2_unit_distance_lab)); 

    return 0;
    }//get_angle_candidates(phi_increment)

int average_candidates(){

    for(int j = 0; j < 6; j++){
        wp_candidate[j] = 0; 
        for(int i = 0; i < number_of_candidates; i++){wp_candidate[j] += vector_wp_lab_candidate[j][i];}
        wp_candidate[j] = wp_candidate[j]/static_cast<float>(number_of_candidates);
        }

    return 0;
    }//average_candidates()

int get_opening_angles(){

    if(number_of_candidates!=0){    
        daw1 = 0; 
        daw2 = 0;
        for(int l = 0; l < number_of_candidates; l++){
            daw1 = daw1 + vector_w1w1[l];
            daw2 = daw2 + vector_w2w2[l];
            }

        daw1 = daw1/static_cast<float>(number_of_candidates);
        daw2 = daw2/static_cast<float>(number_of_candidates);

        //daw1_all.push_back(daw1);
        //daw2_all.push_back(daw2);
        }

    return 0;
    }//get_opening_angles()

int append_values(){

    vector_w1phi_lab_candidate.push_back(w1_lab_candidate.Phi()); 
    vector_w2phi_lab_candidate.push_back(w2_lab_candidate.Phi());

    vector_wp_lab_candidate[0].push_back(w1_lab_candidate.Px()); 
    vector_wp_lab_candidate[1].push_back(w1_lab_candidate.Py()); 
    vector_wp_lab_candidate[2].push_back(w1_lab_candidate.Pz()); 
                       
    vector_wp_lab_candidate[3].push_back(w2_lab_candidate.Px()); 
    vector_wp_lab_candidate[4].push_back(w2_lab_candidate.Py());
    vector_wp_lab_candidate[5].push_back(w2_lab_candidate.Pz()); 

    //vector_wangle_lab[0].push_back(l1w1_angle_lab); 
    //vector_wangle_lab[1].push_back(l2w2_angle_lab); 

    for (int i = 0; i < 6; i++){
        //all_wp_lab_candidate[i].insert(all_wp_lab_candidate[i].end(), vector_wp_lab_candidate[i].begin(), vector_wp_lab_candidate[i].end());
    }
    
    vector_w1w1.push_back(acos(w1_unit_distance_lab.Dot(w1_unit_distance_lab_ref)));
    vector_w2w2.push_back(acos(w2_unit_distance_lab.Dot(w2_unit_distance_lab_ref)));
    vector_w1w2.push_back(acos(w2_unit_distance_lab.Dot(w1_unit_distance_lab_ref)));
    vector_w2w1.push_back(acos(w1_unit_distance_lab.Dot(w2_unit_distance_lab_ref))); 
    

    return 0;
    }//append_values(i,j)

int clear_vectors(){

    vector_w1phi_lab_candidate.clear();
    vector_w2phi_lab_candidate.clear();
    //vector_wangle_lab[0].clear();
    //vector_wangle_lab[1].clear();

    for(int k = 0; k < 6; k++){vector_wp_lab_candidate[k].clear();}

    vector_w1w1.clear();
    vector_w2w2.clear();
    vector_w1w2.clear();
    vector_w2w1.clear();

    return 0;
    }//clear_vectors()

int average_filter(float max_wwe_diff, float angle_limit){
    bool pass_cut = false; // Satisfy both the cuts?
    float w1e_candidate_av; 
    float w2e_candidate_av;
    float cut_angle;
            
    // For events with more than one candidates take the average
    number_of_candidates = vector_w1phi_lab_candidate.size(); // Getting number of events

    while(pass_cut == false){
        //cout << " LOOPING " << endl;
        if (number_of_candidates == 0) {
            pass_cut = true;
            break;
            }
        average_candidates(); // average from the selection that satisfied the cuts
        w1e_candidate_av = sqrt(w_mass*w_mass + wp_candidate[0]*wp_candidate[0] + wp_candidate[1]*wp_candidate[1] + wp_candidate[2]*wp_candidate[2]); // Calculate energy of averaged w1
        w2e_candidate_av = sqrt(w_mass*w_mass + wp_candidate[3]*wp_candidate[3] + wp_candidate[4]*wp_candidate[4] + wp_candidate[5]*wp_candidate[5]); // Calculate energy of averaged w2
        w1_lab_average.SetCoordinates(wp_candidate[0], wp_candidate[1], wp_candidate[2], w1e_candidate_av); // 4-momnetum of averaged w1
        w2_lab_average.SetCoordinates(wp_candidate[3], wp_candidate[4], wp_candidate[5], w2e_candidate_av); // 4-momnetum of averaged w2

        if(abs(w1e_candidate_av + w2e_candidate_av - wwe_lab) < max_wwe_diff){ // Energy cut on averaged results

            //cout << " PASSED ENERGY CUT " << endl;

            l1w1_angle_lab = 0;
            l2w2_angle_lab = 0;
            w1_lab = w1_lab_average;
            w2_lab = w2_lab_average;
            w1e_lab = w1_lab.E(); w1px_lab = w1_lab.Px(); w1py_lab = w1_lab.Py(); w1pz_lab = w1_lab.Pz();
            w2e_lab = w2_lab.E(); w2px_lab = w2_lab.Px(); w2py_lab = w2_lab.Py(); w2pz_lab = w2_lab.Pz();
            w1_beta = sqrt(w1px_lab*w1px_lab + w1py_lab*w1py_lab + w1pz_lab*w1pz_lab)/w1e_lab; 
            w2_beta = sqrt(w2px_lab*w2px_lab + w2py_lab*w2py_lab + w2pz_lab*w2pz_lab)/w2e_lab;
            w1_gamma = 1/sqrt(1-w1_beta*w1_beta);
            w2_gamma = 1/sqrt(1-w2_beta*w2_beta);
            l1w1_angle_lab = acos( 1/w1_beta - w1_gamma * w_mass*(1/w1_beta - w1_beta)/2/l1e_lab_event ); 
            l2w2_angle_lab = acos( 1/w2_beta - w2_gamma * w_mass*(1/w2_beta - w2_beta)/2/l2e_lab_event ); 

            l1w1_angle_lab_candidate = acos(l1_unit_distance_lab.Dot(w1_lab_average.Vect().Unit()));
            l2w2_angle_lab_candidate = acos(l2_unit_distance_lab.Dot(w2_lab_average.Vect().Unit()));

            /*
            cout << "Average l1w1_angle_lab" << l1w1_angle_lab << endl;
            cout << "Average l2w2_angle_lab" << l2w2_angle_lab << endl;
            cout << "Average l1w1_angle_lab_candidate" << l1w1_angle_lab_candidate << endl;
            cout << "Average l2w2_angle_lab_candidate" << l2w2_angle_lab_candidate << endl;
            */

            if( (abs(l1w1_angle_lab_candidate - l1w1_angle_lab) <= angle_limit)&&(abs(l2w2_angle_lab_candidate - l2w2_angle_lab) <= angle_limit) ){

                //cout << " PASSED ANGLE CUT " << endl;

                pass_cut = true;
                
                //vector_w1w1.clear();
                //vector_w2w2.clear();
                //w1_unit_distance_lab = w1_lab_average.Vect().Unit();
                //w2_unit_distance_lab = w2_lab_average.Vect().Unit();
                //vector_w1w1.push_back(acos((w1_unit_distance_lab).Dot(w1_unit_distance_lab_ref)));
                //vector_w2w2.push_back(acos((w2_unit_distance_lab).Dot(w2_unit_distance_lab_ref)));
                
                //number_of_candidates = 1;
                }//if
             
            /*
            l1w1_angle_lab_candidate = acos(l1_unit_distance_lab.Dot(w1_lab_average.Vect().Unit())); 
            l2w2_angle_lab_candidate = acos(l2_unit_distance_lab.Dot(w2_lab_average.Vect().Unit()));
            for(int l = 0; l < number_of_candidates; l++){
                if( (abs(l1w1_angle_lab_candidate - vector_wangle_lab[0][l]) <= angle_limit)&&(abs(l2w2_angle_lab_candidate - vector_wangle_lab[1][l]) <= angle_limit) ){
                    pass_cut = true;
                    vector_w1w1.clear();
                    vector_w2w2.clear();
                    vector_w1w1.push_back(acos((w1_lab_average.Vect().Unit()).Dot(w1_unit_distance_lab_ref)));
                    vector_w2w2.push_back(acos((w2_lab_average.Vect().Unit()).Dot(w2_unit_distance_lab_ref)));
                    }//if
                cut_angle = l; 
                break;
                }//for
                */

            }//if
        if (pass_cut == false){
            float maximum_mag1 = 0;
            float maximum_mag1_value = 0;
            float maximum_mag2 = 0;
            float maximum_mag2_value = 0;
            float maximum_mag1_angle = 0;
            float maximum_mag1_value_angle = 0;
            float maximum_mag2_angle = 0;
            float maximum_mag2_value_angle = 0;
            float w1p_magnitude;
            float w2p_magnitude;
            float w1a_magnitude;
            float w2a_magnitude;
            float maximum_mag_value_angle = 0;
            float maximum_mag = 0;
            for (int c = 0; c < number_of_candidates; c++){ // Selecting largest deviation from average
                //w1p_magnitude = abs(sqrt(w1_lab_average.Px()*w1_lab_average.Px()+w1_lab_average.Py()*w1_lab_average.Py()+w1_lab_average.Pz()*w1_lab_average.Pz())-sqrt(vector_wp_lab_candidate[0][c]*vector_wp_lab_candidate[0][c]+vector_wp_lab_candidate[1][c]*vector_wp_lab_candidate[1][c]+vector_wp_lab_candidate[2][c]*vector_wp_lab_candidate[2][c]));
                //w2p_magnitude = abs(sqrt(w2_lab_average.Px()*w2_lab_average.Px()+w2_lab_average.Py()*w2_lab_average.Py()+w2_lab_average.Pz()*w2_lab_average.Pz())-sqrt(vector_wp_lab_candidate[3][c]*vector_wp_lab_candidate[3][c]+vector_wp_lab_candidate[4][c]*vector_wp_lab_candidate[4][c]+vector_wp_lab_candidate[5][c]*vector_wp_lab_candidate[5][c]));
                w1_angle_diff.SetCoordinates(vector_wp_lab_candidate[0][c],vector_wp_lab_candidate[1][c],vector_wp_lab_candidate[2][c]);
                w2_angle_diff.SetCoordinates(vector_wp_lab_candidate[3][c],vector_wp_lab_candidate[4][c],vector_wp_lab_candidate[5][c]); 
                w1a_magnitude = (w1_angle_diff.Unit()).Dot(w1_lab_average.Vect().Unit());
                w2a_magnitude = (w2_angle_diff.Unit()).Dot(w2_lab_average.Vect().Unit());
                
                /*
                cout << "w1a_magnitude" << w1a_magnitude << endl;
                cout << "w2a_magnitude" << w2a_magnitude << endl;
                */

                // INDEPENDENT CHECKS
                /* 
                if (w1p_magnitude > maximum_mag1_value){
                    maximum_mag1_value = w1p_magnitude;
                    maximum_mag1 = c;
                    }
                if (w2p_magnitude > maximum_mag2_value){
                    maximum_mag2_value = w2p_magnitude;
                    maximum_mag2 = c;
                    }
                if (w1a_magnitude > maximum_mag1_value_angle){
                    maximum_mag1_value_angle = w1a_magnitude;
                    maximum_mag1_angle = c;
                    }
                if (w2a_magnitude > maximum_mag2_value_angle){
                    maximum_mag2_value_angle = w2a_magnitude;
                    maximum_mag2_angle = c;
                    }
                */

                // SIMULTANEOUS CHECKS
                if (w1a_magnitude + w2a_magnitude > maximum_mag_value_angle){
                    maximum_mag_value_angle = w1a_magnitude + w2a_magnitude;
                    maximum_mag = c;
                    }
                /*
                if (w1p_magnitude + w2p_magnitude > maximum_mag_value){
                    maximum_mag_value = w1p_magnitude + w2p_magnitude;
                    maximum_mag = c;
                    }
                */
                }
            // Removing largest deviation from average
            vector_wp_lab_candidate[0].erase(vector_wp_lab_candidate[0].begin()+maximum_mag);
            vector_wp_lab_candidate[1].erase(vector_wp_lab_candidate[1].begin()+maximum_mag);
            vector_wp_lab_candidate[2].erase(vector_wp_lab_candidate[2].begin()+maximum_mag);
            vector_wp_lab_candidate[3].erase(vector_wp_lab_candidate[3].begin()+maximum_mag);
            vector_wp_lab_candidate[4].erase(vector_wp_lab_candidate[4].begin()+maximum_mag);
            vector_wp_lab_candidate[5].erase(vector_wp_lab_candidate[5].begin()+maximum_mag);
            vector_w1w1.erase(vector_w1w1.begin()+maximum_mag);
            vector_w2w2.erase(vector_w2w2.begin()+maximum_mag);
            //vector_wangle_lab[0].erase(vector_wangle_lab[0].begin()+maximum_mag1);
            //vector_wangle_lab[1].erase(vector_wangle_lab[1].begin()+maximum_mag2);
            number_of_candidates -= 1;
            }

        }
    return 0;
    }

//***************************************************// MAIN //***************************************************//
string run_simulation(const int event_count, float max_wwe_diff, float angle_limit, const int max_phi, float smearing_extent, int i, float parameter, float parameter_count) {

    // File directory
    TFile * f = TFile::Open("/Users/swathi/Desktop/Final_GP_S/Simulation_Output.root","RECREATE");
    TTree * t = new TTree("w1w2","W1W2");

    int success_event_count = 0; 
    int non_success_event_count = 0;
    successRate = 0;

    load_data();
    create_directory(t);

	for(int event = 0; event < event_count; event++){ 
        reference_variables(event); 
        fill_lab_values();
	    ww_boost_to_com(smearing_extent);

		for( int i = 0; i < wwe_range; i++ ){

			for( size_t j = 0; j < wwp_range; j++ ){

                wwe_com = w_mass*2 + i*2*(wwe_lab-w_mass)/(wwe_range-1); // WW Energy incrementally increasing
                wwp_com = -sqrt(wwe_com*wwe_com/4 - w_mass*w_mass) + 2*j*sqrt(wwe_com*wwe_com/4 - w_mass*w_mass)/(wwp_range-1); // WW Z-Momentum incrementally increasing
                w1w2_boost_to_lab();
                wwe_diff = abs(w1e_lab + w2e_lab - wwe_lab);

				if(wwe_diff < max_wwe_diff){ // Energy cut
                    
					for( int k = 0; k < max_phi; k++ ){

                        phi_increment = 2*TMath::Pi()*(k+1)/max_phi - TMath::Pi();
                        get_angle_candidates(phi_increment);

                        if( (abs(l1w1_angle_lab_candidate - l1w1_angle_lab) <= angle_limit)&&(abs(l2w2_angle_lab_candidate - l2w2_angle_lab) <= angle_limit) ){ // Angle cut

                            /*
                            cout << "Single l1w1_angle_lab" << l1w1_angle_lab << endl;
                            cout << "Single l2w2_angle_lab" << l2w2_angle_lab << endl;
                            cout << "Single l1w1_angle_lab_candidate" << l1w1_angle_lab_candidate << endl;
                            cout << "Single l2w2_angle_lab_candidate" << l2w2_angle_lab_candidate << endl;
                            */

                            append_values();

                            } // Angle cut
                        }//k
				    }// Energy cut
			    }//j
            }//i

        //average_candidates();
        //get_opening_angles();
        average_filter(max_wwe_diff, angle_limit);
        get_opening_angles();

        // Counting number of successful candidates
		if(number_of_candidates == 0){
            non_success_event_count += 1;
            } 
        else{
            success_event_count += 1;
            t->Fill(); // Outputing data into tree
            }

        cout << "Event " << event + 1 << " complete. Number of candidates: " << number_of_candidates << endl;

        //t->Fill(); // Outputing data into tree

        clear_vectors();

	    }// Events
    
    t->Write();

	//Success rate of the reconstruction
	successRate = ( (float)success_event_count / (float)event_count ) * 100.;
	cout << "***************************************************************************" << endl; 
	cout << "Reconstruction success: " << successRate << "%" << endl;
    cout << "***************************************************************************" << endl; 

    max_events = event_count/2; // Turn in to maximum value of plot + padding

    successRate_out.push_back( successRate );
    //daw1_all_out.push_back( std::accumulate(daw1_all.begin(), daw1_all.end(), 0.0) / daw1_all.size() );
    //daw2_all_out.push_back( std::accumulate(daw2_all.begin(), daw2_all.end(), 0.0) / daw2_all.size() );
    //daw1_all.clear();
    //daw2_all.clear();
    
    string print_location = format_plots(max_events);
    //scatter_plot(print_location,all_wp_lab_candidate,wp_ref);

    return print_location;
    }//run_simulation(event_count,max_wwe_diff,angle_limit,max_phi,smearing_extent)


int loop_simulation(float parameter_start, float parameter_end, int parameter_count){

    // Default Modifiables
    const int event_count = 500;
    float max_wwe_diff = 10; 
    float angle_limit = 0.001;
    int max_phi = 360; 
    float smearing_extent = 0;

    for( int i = 0; i < parameter_count; i++ ){

        float parameter_increase = (parameter_end-parameter_start)/parameter_count;
        float parameter = parameter_start + i*parameter_increase;

        max_wwe_diff = parameter; // ASSIGN WITH PARAMETER TO MODIFY

        auto start = std::chrono::high_resolution_clock::now();
        print_location = run_simulation(event_count,max_wwe_diff,angle_limit,max_phi,smearing_extent,i,parameter,parameter_count); 
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

        ofstream plotfile(print_location+"plot_info.txt");

        plotfile << "Plot information for file: " << print_location << "\n\n";

        plotfile << "--- Trial information ---\n";
        plotfile << "Trial number: " << i+1 << " out of " << parameter_count << "\n";
        plotfile << "Elapsed time: " << duration.count() << " ms" << "\n";
        plotfile << "Parameter value: " << parameter << "\n";
        plotfile << "Success rate: " << successRate << " %\n\n";

        plotfile << "--- Parameter initialisation ---\n";
        plotfile << "Number of simulated events: " << event_count << "\n";
        plotfile << "Energy difference upper cut limit: " << max_wwe_diff << "\n";
        plotfile << "Angle cut limit: " << angle_limit << "\n";
        plotfile << "Azimuthal angle range: " << max_phi << "\n";
        plotfile << "Smearing extent: " << smearing_extent << "\n\n";

        repeats.push_back(parameter);
        time_out.push_back(duration.count());

        }

    //print(daw1_all_out);
    //print(daw2_all_out);
    print(time_out);
    print(successRate_out);
    print(repeats);
    
    format_plots_averages(print_location,repeats,successRate_out,time_out);

    cout << "***************************************************************************" << endl; 
	cout << " TASK COMPLETE " << endl;
    cout << "***************************************************************************" << endl; 

    return 0;
    }//loop_simulation(parameter_start,parameter_end,parameter_count)
