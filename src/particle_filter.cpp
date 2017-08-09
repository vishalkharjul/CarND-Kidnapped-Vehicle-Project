/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"

using namespace std;



void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
    num_particles = 100;
    default_random_engine gen;
    double std_x, std_y, std_psi; // Standard deviations for x, y, and psi

    // TODO: Set standard deviations for x, y, and psi.
    std_x=std[0];
    std_y=std[1];
    std_psi=std[2];


    // This line creates a normal (Gaussian) distribution for x.
    normal_distribution<double> dist_x(x, std_x);

    // TODO: Create normal distributions for y and psi.
    normal_distribution<double> dist_y(y, std_y);
    normal_distribution<double> dist_psi(theta, std_psi);

    for (int i = 0; i < num_particles; ++i) {
        //double sample_x, sample_y, sample_psi;
        Particle sample;
        // TODO: Sample  and from these normal distrubtions like this:
        sample.id=i;
        sample.x = dist_x(gen);
        sample.y = dist_y(gen);
        sample.theta = dist_psi(gen);
        sample.weight = 1.0;
        weights.push_back(1.0);
        particles.push_back(sample);

        // where "gen" is the random engine initialized earlier (line 18).

        // Print your samples to the terminal.
        // cout << "Sample " << i + 1 << " " << sample_x << " " << sample_y << " " << sample_psi << endl;
    }
    is_initialized = true;



}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
    double std_x, std_y, std_psi; // Standard deviations for x, y, and psi

    std_x=std_pos[0];
    std_y=std_pos[1];
    std_psi=std_pos[2];

    default_random_engine gen;
    normal_distribution<double> dist_x(0, std_x);
    normal_distribution<double> dist_y(0, std_y);
    normal_distribution<double> dist_psi(0, std_psi);

    for (auto& particle:particles) {
        //double sample_x, sample_y, sample_psi;

        if(fabs(yaw_rate)>0.001) {
            double vyaw_rate = velocity/yaw_rate;
            particle.x = particle.x + vyaw_rate * (sin(particle.theta + (yaw_rate*delta_t))-sin(particle.theta));
            particle.y = particle.y + vyaw_rate * (cos(particle.theta)-cos(particle.theta+(yaw_rate*delta_t)));
            particle.theta = particle.theta + (delta_t*yaw_rate);

        } else{
            particle.x += velocity * delta_t * cos(particle.theta);
            particle.y += velocity * delta_t * sin(particle.theta);
        }


        particle.x = particle.x + dist_x(gen);
        particle.y = particle.y + dist_y(gen);
        particle.theta = particle.theta+dist_psi(gen);
        // where "gen" is the random engine initialized earlier (line 18).

    }


}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

}

double ParticleFilter::weight_prob(double x, double y, float xm, float ym, double std_landmark[]) {
    double xpart = -0.5*(x - xm)*(x - xm) / (std_landmark[0] * std_landmark[0]);
    double ypart = -0.5*(y - ym)*(y - ym) / (std_landmark[1] * std_landmark[1]);
    return exp(xpart + ypart) / (2 * M_PI*std_landmark[0] * std_landmark[1]);
}



void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		std::vector<LandmarkObs> observations, Map map_landmarks) {
    // TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
    //   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
    // NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
    //   according to the MAP'S coordinate system. You will need to transform between the two systems.
    //   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
    //   The following is a good resource for the theory:
    //   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
    //   and the following is a good resource for the actual equation to implement (look at equation
    //   3.33
    //   http://planning.cs.uiuc.edu/node99.html

    weights.clear();
    Map selected_map;
    selected_map.landmark_list.clear();

    for (auto& particle : particles) // loop every particle
    {
        double par_x = particle.x;
        double par_y = particle.y;
        double par_theta = particle.theta;


        for (int j = 0; j < map_landmarks.landmark_list.size(); ++j)  // loop for very landmarks
        {
            if (dist(par_x,par_y,map_landmarks.landmark_list[j].x_f,map_landmarks.landmark_list[j].y_f)<sensor_range)

            {
                selected_map.landmark_list.push_back(map_landmarks.landmark_list[j]);
            }
        }

        // transform coordinate system
        particle.weight = 1.0;
        LandmarkObs tmp_ob;
        for (auto& tmp_ob : observations) // loop for very observations
        {
            double temp_x = tmp_ob.x * cos(par_theta) - tmp_ob.y * sin(par_theta) + par_x;
            double temp_y = tmp_ob.x * sin(par_theta) + tmp_ob.y * cos(par_theta) + par_y;

            double min = 10000.0;
            int min_id = 0;

            // association
            double tmp_dist;
            tmp_ob.id = 0;
            for (int k = 0; k < selected_map.landmark_list.size(); ++k)  // loop for every possible landmarks
            {
                tmp_dist = dist(selected_map.landmark_list[k].x_f ,selected_map.landmark_list[k].y_f, temp_x,temp_y);
                if (tmp_dist < min)
                {
                    min = tmp_dist;
                    //tmp_ob.id = selected_map.landmark_list[k].id_i;
                    min_id = selected_map.landmark_list[k].id_i;
                }
            }

            if (observations.size() == 0) {
                particle.weight = 0.0;
            }
            else
            {
                particle.weight *= weight_prob(temp_x, temp_y, map_landmarks.landmark_list[min_id-1].x_f, map_landmarks.landmark_list[min_id-1].y_f, std_landmark);
            }


        }
        weights.push_back(particle.weight);
        selected_map.landmark_list.clear();
    }

}


void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
    std::vector<Particle> selected_particles;
    std::default_random_engine generator;

    std::discrete_distribution<> distribution (weights.begin(), weights.end());

    for (int i = 0; i < num_particles; ++i)
    {
        int idx = distribution(generator);
        selected_particles.push_back(particles[idx]);
    }


    particles = selected_particles;

}

Particle ParticleFilter::SetAssociations(Particle particle, std::vector<int> associations, std::vector<double> sense_x, std::vector<double> sense_y)
{
	//particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
	// associations: The landmark id that goes along with each listed association
	// sense_x: the associations x mapping already converted to world coordinates
	// sense_y: the associations y mapping already converted to world coordinates

	//Clear the previous associations
	particle.associations.clear();
	particle.sense_x.clear();
	particle.sense_y.clear();

	particle.associations= associations;
 	particle.sense_x = sense_x;
 	particle.sense_y = sense_y;

 	return particle;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
