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
#include <map>

#include "particle_filter.h"

using namespace std;

//#define FILEWRITE
void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of
	//   x, y, theta and their uncertainties from GPS) and all weights to 1.
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
#ifdef FILEWRITE
	ofstream file;
	file.open ("codebind.txt");
	file << "Entering Init \n";
	file.close();
#endif
	num_particles = 10; // Less than that causes wrong estimations, more than that is too much computation
	double std_x = std[0];
	double std_y = std[1];
	double std_theta = std[2];

	std::default_random_engine generator;

	std::normal_distribution<double> dist_x(x,std_x);
	std::normal_distribution<double> dist_y(y,std_y);
	std::normal_distribution<double> dist_theta(theta,std_theta);

	weights.resize(num_particles,1.0);
	Particle p;
	particles.resize(num_particles,p);
	if (!is_initialized)
	{
		for (int i = 0;i < num_particles;i++)
		{
			particles[i].x = dist_x(generator);
			particles[i].y = dist_y(generator);
			particles[i].theta = dist_theta(generator);
			particles[i].weight = 1;
			particles[i].id = i + 1;
			weights[i] = 1;
#ifdef FILEWRITE
			file.open ("codebind.txt",std::fstream::app);
			file<<"Particle No. "<<i + 1<<"  "<<"x= "<<particles[i].x<<"	"<<"y= "<<particles[i].y<<"	"<<"Theta= "<<particles[i].theta<<"	"<<"Weight= "<<particles[i].weight <<endl;
			file.close();
#endif
		}
		is_initialized = true;

	}
#ifdef FILEWRITE
	file.open ("codebind.txt",std::fstream::app);
	file << "Exiting Init \n";
	file.close();
#endif

}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
#ifdef FILEWRITE
	ofstream file;
	file.open ("codebind.txt",std::fstream::app);
	file << "Entering Prediction \n";
	file.close();
#endif
	double std_x = std_pos[0];
    double std_y = std_pos[1];
    double std_theta = std_pos[2];

    std::default_random_engine generator2;

    std::normal_distribution<double> dist_x(0,std_x);
    std::normal_distribution<double> dist_y(0,std_y);
    std::normal_distribution<double> dist_theta(0,std_theta);

    for (int i = 0;i < num_particles; i++)
    {
    	double x = particles[i].x;
    	double y = particles[i].y;
    	double theta = particles[i].theta;
    	double x_new;
    	double y_new;
    	double theta_new;

    	if (abs(yaw_rate) < 0.001)
    	{
    		x_new = x + velocity*cos(theta)*delta_t;
    		y_new = y + velocity*sin(theta)*delta_t;
    		theta_new = theta;
    	}
    	else
    	{
    		x_new = x + velocity/yaw_rate * (sin(theta + yaw_rate*delta_t) - sin(theta));
    		y_new = y + velocity/yaw_rate * (cos(theta) - cos(theta+yaw_rate*delta_t));
    		theta_new = theta + yaw_rate*delta_t;
    	}
    	particles[i].x = x_new + dist_x(generator2);
    	particles[i].y = y_new + dist_y(generator2);
    	particles[i].theta = theta_new + dist_theta(generator2);

#ifdef FILEWRITE
    	file.open ("codebind.txt",std::fstream::app);
    	file<<"Particle No. "<<i + 1<<"  "<<"x= "<<particles[i].x<<"	"<<"y= "<<particles[i].y<<"	"<<"Theta= "<<particles[i].theta<<"	"<<"Weight= "<<particles[i].weight <<endl;
    	file.close();
#endif
    }
#ifdef FILEWRITE
    	  		file.open ("codebind.txt",std::fstream::app);
    	  		file << "Exiting prediction \n";
    	  		file.close();
#endif
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to
	//   implement this method and use it as a helper during the updateWeights phase.
#ifdef FILEWRITE
	ofstream file;
	file.open ("codebind.txt",std::fstream::app);
	file << "Entering data Association \n";
	file.close();
#endif
	double nearest_dist,o_x,o_y,pred_x,pred_y;
	for (int i = 0;i < observations.size();i++)
	{
		nearest_dist = std::numeric_limits<double>::max();
		o_x = observations[i].x;
		o_y = observations[i].y;
		for (int j = 0; j < predicted.size();j++)
		{
			pred_x = predicted[j].x;
			pred_y = predicted[j].y;
			if (dist(o_x,o_y,pred_x,pred_y) < nearest_dist)
			{
				nearest_dist = dist(o_x,o_y,pred_x,pred_y);
				observations[i].id = predicted[j].id;
			}
		}
	}
#ifdef FILEWRITE
		  		file.open ("codebind.txt",std::fstream::app);
		  		file << "Exiting dataAssociation \n";
		  		file.close();
#endif

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

#ifdef FILEWRITE
	ofstream file;
	file.open ("codebind.txt",std::fstream::app);
	file << "Entering updateWeights \n";
	file.close();
#endif
	double std_x = std_landmark[0];
	double std_y = std_landmark[1];

	std::vector<LandmarkObs> global_observations;
	std::vector<LandmarkObs> inrange_landmarks;

	weights.clear();
	for (int i = 0; i < num_particles; i++)
	{
		double p_x = particles[i].x;
		double p_y = particles[i].y;
		double p_theta = particles[i].theta;
		double global_x, global_y,obs_dist,landmark_dist,l_x,l_y;
		int l_id;

		global_observations.clear();
		LandmarkObs g_obs,g_landmark;

		//double weight = 1;
		particles[i].weight = 1; // Setting weights bck to 1

		for (int j = 0; j < observations.size(); j++)
		{
			global_x = (observations[j].x*cos(p_theta) - observations[j].y*sin(p_theta)) + p_x; // Transformation into MAP coordinate system
			global_y = (observations[j].x*sin(p_theta) + observations[j].y*cos(p_theta)) + p_y;
			//obs_dist = dist(0,0,global_x,global_y);
			if (true) // obs_dist < sensor_range
			{
				g_obs.x = global_x;
				g_obs.y = global_y;
				g_obs.id = observations[j].id;
				global_observations.push_back(g_obs);
			}
		}
		inrange_landmarks.clear();

		for (int k = 0; k < map_landmarks.landmark_list.size();k++)
		{
			l_x = map_landmarks.landmark_list[k].x_f;
			l_y = map_landmarks.landmark_list[k].y_f;
			landmark_dist = dist(p_x,p_y,l_x,l_y);
			if (landmark_dist < sensor_range) // Considering only the landmarks with in sensor range
			{
				g_landmark.x = l_x;
				g_landmark.y = l_y;
				g_landmark.id = map_landmarks.landmark_list[k].id_i;
				inrange_landmarks.push_back(g_landmark);
			}
		}

		dataAssociation(inrange_landmarks,global_observations); // Associate ovservations to landmarks

		double o_x, o_y, landmark_x, landmark_y,dx,dy;
		int o_id, landmark_id;

		for (int m = 0; m < global_observations.size(); m++)
		{
			o_x = global_observations[m].x;
			o_y = global_observations[m].y;
			o_id =global_observations[m].id;

			for (int n = 0; n < inrange_landmarks.size(); n++)
			{
				if (inrange_landmarks[n].id == o_id)
				{
					landmark_x = inrange_landmarks[n].x;
					landmark_y = inrange_landmarks[n].y;
#ifdef FILEWRITE
					file.open ("codebind.txt",std::fstream::app);
					file<<"Landmark ID. "<<inrange_landmarks[n].id <<endl;
					file<<"Transformed Obs. x= "<<o_x<<  "	Transformed Obs. y= "<<o_y<<endl;
					file<<"Landmark x= "<<landmark_x<<  "	Landmark y= "<<landmark_y<<endl;
					file.close();
#endif
					break; // We have found the associated landmark, no need to loop further
				}

			}
			dx = o_x - landmark_x;
			dy = o_y - landmark_y;
#ifdef FILEWRITE
			file.open ("codebind.txt",std::fstream::app);
			file<<"Difference x = "<<dx<<  "	Difference y = "<<dy<<endl;
#endif
			double d1 = (dx*dx)/(2.0*std_x*std_x);
			double d2 = (dy*dy)/(2.0*std_y*std_y);
			double p = (0.5/(M_PI*std_x*std_y))*exp(-(d1 + d2)); // Bivariate pdf of an observation-landmark pair 
#ifdef FILEWRITE
			file<<"Instantaneous Probability = "<<p<< endl;
#endif
			particles[i].weight = particles[i].weight * p; // Cummulative weight
			//weight = weight*p;
#ifdef FILEWRITE
			file<<"Cumulative Probability = "<<particles[i].weight << endl;
			file.close();
#endif
		}
		//particles[i].weight = particles[i].weight*weight;
		weights.push_back(particles[i].weight);
#ifdef FILEWRITE
		file.open ("codebind.txt",std::fstream::app);
		file<<"Final probability = "<<particles[i].weight<<endl;
		file.close();
		file.open ("codebind.txt",std::fstream::app);
		file<<"Particle No. "<<i + 1<<"  "<<"x= "<<particles[i].x<<"	"<<"y= "<<particles[i].y<<"	"<<"Theta= "<<particles[i].theta<<"	"<<"Weight= "<<particles[i].weight <<endl;
		file.close();
#endif
	}
#ifdef FILEWRITE
		  		file.open ("codebind.txt",std::fstream::app);
		  		file << "Exiting updateWeights \n";
		  		file.close();
#endif
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight.
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
	//std::discrete_distribution<int> dist_ind(particles[0].weight, particles[num_particles -1].weight);
#ifdef FILEWRITE
	ofstream file;
	file.open ("codebind.txt",std::fstream::app);
	file << "Entering Resample \n";
	file.close();
#endif
	std::discrete_distribution<int> dist_ind(weights.begin(), weights.end());
	std::vector<Particle> new_particles;
	std::default_random_engine generator3;

	new_particles.clear();
	int index;
	for (int i = 0; i < num_particles; i++)
	{
		index = dist_ind(generator3); // Indices with higher weights should be choosen more often
#ifdef FILEWRITE
		file.open ("codebind.txt",std::fstream::app);
		file<<"Particle Id.  = "<<i + 1<<"  Random chosen index. "<<index<<endl;
		file.close();
#endif
		Particle p;
		p = particles[index];
		p.id = i + 1;
		new_particles.push_back(p);
	}
	particles = new_particles; // updated particle list
#ifdef FILEWRITE
		  		file.open ("codebind.txt",std::fstream::app);
		  		file << "Exiting Resample \n";
		  		file << "Cycle completed \n";
		  		file << "\n";
		  		file.close();
#endif
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

