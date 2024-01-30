#include <iostream>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <time.h>
#include <omp.h>
// function to initialise the grid
void initialiseGrid(std::vector<std::vector<int>>& grid, double p, int N){
    
    for(auto& row : grid) {
        for (int& cell : row){
            double randValue = static_cast<double>(rand())/ RAND_MAX;
            if (randValue < p) {
                cell = 1;
                }
            else {
                cell = 0;
            }
        }
    }

    int topBurning = 0;
    for (int i = 0; i<N ; ++i) {
    
        // set the whole top row alight
        if (grid[0][i] == 1){
            grid[0][i] = 2;
            

        }
            
    }
    
    
}


struct SimulationResult {
    int steps;
    int fireReachedBottom;
    double timeElapsed;
};

//  function to simulate the forest fire
SimulationResult simulateForestFire(std::vector<std::vector<int>>& grid, int N, int& fireReachedBottom, int& steps) {
    steps = 0;
    bool fireStillBurning = true;
    // start timer - only timing the main loop
    struct timespec start, finish;
    clock_gettime(CLOCK_MONOTONIC, &start);
    while (fireStillBurning == true) {
        // Create a copy of the grid to store the next state
        std::vector<std::vector<int>> nextGrid = grid;
        // 2 is burning, 1 is alive, 0 is dead
        #pragma omp parallel for schedule(dynamic)
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                if(grid[N-1][j] ==2){
                    fireReachedBottom = 1;
                }
                if (grid[i][j] == 1) {  
                    int neighborsX[] = {0, 0, -1, 1};
                    int neighborsY[] = {-1, 1, 0, 0};

                    for (int k = 0; k < 4; ++k) {
                        int ni = i + neighborsX[k];
                        int nj = j + neighborsY[k];
                        // looping over the grid to update trees that are burning
                        if (ni >= 0 && ni < N && nj >= 0 && nj < N) {
                            if (grid[ni][nj] == 2 && grid[i][j] == 1) {
                                nextGrid[i][j] = 2;
                            }

                        }
                    }
                }
                // the trees that were burning are now burnt
                if (grid[i][j] == 2) {
                    nextGrid[i][j] = 0;
                }
                // trees that were alive, and not defined as burning, remain alive           
                if (grid[i][j] == 1 && nextGrid[i][j] != 2) {
                    nextGrid[i][j] = 1;
                }
                
            }
        }
        
        // updating the grid to equal the new grid
        grid = nextGrid;

        // Checking if the system is still dynamic
        steps++;

        
        fireStillBurning = false;  // Reset the flag before checking the grid
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                if (grid[i][j] == 2) {
                    fireStillBurning = true;
                    
                    break;
                }
            }
            
        }
        
        if (!fireStillBurning) {
            std::cout << "burnt out";

            break;
        }

    }

    // end timer
    clock_gettime(CLOCK_MONOTONIC, &finish);
    double timeElapsed = (finish.tv_sec - start.tv_sec);
    timeElapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
    
    return {steps, fireReachedBottom, timeElapsed};

}

int main() {
    int N = 500; // Size of grid
    
    int M = 10; // Number of simulations
    double stepSize = 0.1;

    // Output file for CSV results
    std::ofstream csvFile("results.csv");
    csvFile << "p, Average Steps, Average Time Elapsed, Average Fire Reached Bottom" << std::endl;

    std::vector<SimulationResult> results;
    int fireReachedBottom = 0;
    int steps = 0;
    
    for (double p = 0.0; p <= 1.0; p += stepSize) {

        for (int m = 0; m < M; ++m) {
            std::vector<std::vector<int>> grid(N, std::vector<int>(N, 0));
            initialiseGrid(grid, p, N);
            fireReachedBottom = 0;
            
            SimulationResult result = simulateForestFire(grid, N, fireReachedBottom, steps);
            #pragma omp critical
            results.push_back(result);

            // Output time elapsed and steps for each simulation
            std::cout << "Simulation " << m + 1 << ": Time elapsed = " << result.timeElapsed << " s, Steps = " << result.steps << std::endl;
        }
    
        // Calculate averages
        double averageSteps = 0.0;
        double averageTimeElapsed = 0.0;
        double averageFireReachedBottom = 0.0;

        for (const auto& result : results) {
            averageSteps += result.steps;
            averageTimeElapsed += result.timeElapsed;
            averageFireReachedBottom += result.fireReachedBottom;
        }

        averageSteps /= M;
        averageTimeElapsed /= M;
        averageFireReachedBottom /= M;
        
       
        // Output averages to console
        std::cout << "For p = " << p << ": Average steps = " << averageSteps << ", Average time elapsed = " << averageTimeElapsed << " s, Average fire reached bottom = " << averageFireReachedBottom << std::endl;

        // Write results to CSV file
        csvFile << p << ", " << averageSteps << ", " << averageTimeElapsed << ", " << averageFireReachedBottom << std::endl;

        // Clear the results vector for the next iteration
        results.clear();
    }

    csvFile.close();
    return 0;
}



