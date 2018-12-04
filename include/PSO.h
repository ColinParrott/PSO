#pragma once
#ifndef PSO_PSO_H
#define PSO_PSO_H

#include <vector>
#include <random>
#include <functional>

static std::random_device rd;   // random device engine, usually based on /dev/random on UNIX-like systems
static std::mt19937_64 rng(rd());  // initialize Mersennes' twister using rd to generate the seed

class PSO {
private:
    struct Particle {
        std::vector<double> velocity;
        std::vector<double> position;
        std::vector<double> best_position;
        double fitness = INT_MAX;
        double best_fitness = INT_MAX;
    };

    struct Swarm {
        double global_best_fitness;
        std::vector<double> global_best_position;
        std::vector<Particle> particles;
    };


    double omega = 0.7;
    double alpha1 = 2;
    double alpha2 = 2;
    Swarm swarm;
    std::vector<double> max_velocity;


    const int MAX_ITERATIONS;
    const int DIMENSIONS;
    const int NUM_PARTICLES;
    const std::vector<std::vector<double>> BOUNDS;
    std::function<double(std::vector<double>)> benchmarkFunction;

    void initSwarm();
    double generateRandomDouble(double min, double max);
    void updateParticle(Particle &particle);
    std::vector<double> sumVectors(const std::vector<double> &vector1, const std::vector<double> &vector2, const std::vector<double> &vector3);
    void elementWiseMultiply(std::vector<double> &vec, std::vector<double> &r);
    std::vector<double> getRandomNoise(double a, int length);
    double calculateFitness(Particle &particle);
    void updatePosition(Particle &particle);

public:
    PSO(std::function<double(std::vector<double>)> function, std::vector<std::vector<double>> bounds, int dimensions, int num_particles, int max_iter);
    std::vector<double> run();
    void printPopulation();


};


#endif //PSO_PSO_H
