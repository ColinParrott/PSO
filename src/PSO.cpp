#include <utility>
#include <iostream>
#include "../include/PSO.h"


using namespace std;

PSO::PSO(function<double(vector<double>)> function, vector<vector<double>> bounds, int dimensions,
         int num_particles, int max_iter)
        : DIMENSIONS(dimensions), MAX_ITERATIONS(max_iter), NUM_PARTICLES(num_particles), BOUNDS(move(bounds)) {
    benchmarkFunction = move(function);

}


vector<double> PSO::run() {

    initSwarm();
    for (int gen = 0; gen < MAX_ITERATIONS; gen++) {
        for (int i = 0; i < swarm.particles.size(); i++) {
            updateParticle(swarm.particles[i]);
        }
    }

    return swarm.global_best_position;
}

void PSO::initSwarm() {

    swarm = Swarm();
    for (int i = 0; i < DIMENSIONS; i++) {
        swarm.global_best_position.push_back(0.0);

        double range = BOUNDS[i][1] - BOUNDS[i][0];
        max_velocity.push_back(0.1 * range);
    }

    for (int i = 0; i < NUM_PARTICLES; i++) {
        Particle p;
        p.velocity = vector<double>(DIMENSIONS);
        p.fitness = INT_MAX;

        for (int j = 0; j < DIMENSIONS; j++) {
            double min = BOUNDS[j][0];
            double max = BOUNDS[j][1];
            p.position.push_back(generateRandomDouble(min, max));
        }

        p.best_position = p.position;
        swarm.particles.push_back(p);
    }
    
}

double PSO::generateRandomDouble(double min, double max) {
    uniform_real_distribution<double> dist(min, max);
    return dist(rng);
}

void PSO::updateParticle(Particle &particle) {
    vector<double> a1r1 = getRandomNoise(alpha1, particle.velocity.size());
    vector<double> a2r2 = getRandomNoise(alpha2, particle.velocity.size());

    // omega & velocity
    vector<double> momentum;
    momentum.reserve(particle.velocity.size());
    for (int i = 0; i < particle.velocity.size(); i++) {
        momentum.push_back(omega * particle.velocity[i]);
    }

    // particle best pos - current pos
    vector<double> particleBestDifference;
    particleBestDifference.reserve(particle.position.size());
    for (int i = 0; i < particle.position.size(); i++) {
        particleBestDifference.push_back(particle.best_position[i] - particle.position[i]);
    }

    // global best pos - current pos
    vector<double> globalBestDifference;
    globalBestDifference.reserve(particle.position.size());
    for (int i = 0; i < particle.position.size(); i++) {
        globalBestDifference.push_back(swarm.global_best_position[i] - particle.position[i]);
    }

    elementWiseMultiply(particleBestDifference, a1r1);
    elementWiseMultiply(globalBestDifference, a2r2);

    particle.velocity = min(max_velocity, sumVectors(momentum, particleBestDifference, globalBestDifference));
    updatePosition(particle);
    calculateFitness(particle);
}

std::vector<double> PSO::sumVectors(const std::vector<double> &vector1, const std::vector<double> &vector2,
                                    const std::vector<double> &vector3) {
    vector<double> result;
    result.reserve(vector1.size());
    for (int i = 0; i < vector1.size(); i++) {
        result.push_back(vector1[i] + vector2[i] + vector3[i]);
    }

    return result;
}

void PSO::elementWiseMultiply(std::vector<double> &vec, std::vector<double> &r) {
    for (int i = 0; i < vec.size(); i++) {
        vec[i] = vec[i] * r[i];
    }
}

std::vector<double> PSO::getRandomNoise(double a, int length) {
    vector<double> r;
    r.reserve(length);
    for (int i = 0; i < length; i++) {
        r.push_back(a * generateRandomDouble(0.0f, 1.0f));
    }

    return r;
}

double PSO::calculateFitness(PSO::Particle &particle) {
    particle.fitness = benchmarkFunction(particle.position);

    if (particle.fitness < particle.best_fitness) {
        particle.best_position = particle.position;
        particle.best_fitness = particle.fitness;
    }

    if (particle.best_fitness < swarm.global_best_fitness) {
        swarm.global_best_fitness = particle.best_fitness;
        swarm.global_best_position = particle.best_position;
    }


    return particle.fitness;
}

void PSO::updatePosition(PSO::Particle &p) {
    for (int i = 0; i < p.position.size(); i++) {
        p.position[i] = p.position[i] + p.velocity[i];
        if (p.position[i] > BOUNDS[i][1]) {
            p.position[i] = BOUNDS[i][1];
        } else if (p.position[i] < BOUNDS[i][0]) {
            p.position[i] = BOUNDS[i][0];
        }
    }
}

void PSO::printPopulation() {
    cout << "hi" << endl;
    for(Particle p  : swarm.particles){
        for(double d : p.position){
            cout << d << ' ';
        }
        cout << endl;
    }
}
