#include <iostream>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <time.h>
#include <omp.h>

// function to initialise the grid
void initialiseGrid(std::vector<std::vector<int>>& grid, double p, int N) {
    for (auto& row : grid) {
        for (int& cell : row) {
            double randValue = static_cast<double>(rand()) / RAND_MAX;
            if (randValue < p) {
                cell = 1;
            } else {
                cell = 0;
            }
        }
    }


    int topBurning = 0;
    for (int i = 0; i < N; ++i) {
        // set the whole top row alight
        if (grid[0][i] == 1) {
            grid[0][i] = 2;
        }
    }
}
// funciton for animation
void printGridToFile(const std::vector<std::vector<int>>& grid, std::ofstream& outputFile) {
    for (const auto& row : grid) {
        for (int cell : row) {
            outputFile << cell << ' ';
        }
        outputFile << '\n';
    }
    outputFile << '\n';
}

struct SimulationResult {
    int steps;
    int fireReachedBottom;
    double timeElapsed;
};
//  function to simulate the forest fire
SimulationResult simulateForestFire(std::vector<std::vector<int>>& grid, int N, int& fireReachedBottom, int& steps, int windDirection) {
    steps = 0;
    bool fireStillBurning = true;


    // start timer - only timing the main loop
    struct timespec start, finish;
    clock_gettime(CLOCK_MONOTONIC, &start);

    while (fireStillBurning == true) {
        // Create a copy of the grid to store the next state
        std::vector<std::vector<int>> nextGrid = grid;

        // 2 is burning, 1 is alive, 0 is dead
        // Inside the loop that updates the grid
        #pragma omp parallel for schedule(dynamic)
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                if (grid[N - 1][j] == 2) {
                    fireReachedBottom = 1;
                }
                if (grid[i][j] == 2) {
                    int neighborsX[] = {0, 0, -1, 1};
                    int neighborsY[] = {-1, 1, 0, 0};

                    for (int k = 0; k < 4; ++k) {
                        int ni = i + neighborsX[k];
                        int nj = j + neighborsY[k];

                        // Check if wind affects the fire
                        if (windDirection == k && ni >= 0 && ni < N && nj >= 0 && nj < N) {
                            int nni = ni + neighborsX[k];
                            int nnj = nj + neighborsY[k];

                            if (nni >= 0 && nni < N && nnj >= 0 && nnj < N) {
                                // If there are trees in the cells, set them alight
                                if (grid[ni][nj] == 1) {
                                    nextGrid[ni][nj] = 2;
                                }
                                if (grid[nni][nnj] == 1) {
                                    nextGrid[nni][nnj] = 2;
                                }
                            }
                        } else {
                            // looping over the grid to update trees that are burning
                            if (ni >= 0 && ni < N && nj >= 0 && nj < N) {
                                if (grid[ni][nj] == 1) {
                                    nextGrid[ni][nj] = 2;
                                }
                            }
                        }
                    }
                }
            }
        }

        // the trees that were burning are now burnt
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                if (grid[i][j] == 2) {
                    nextGrid[i][j] = 0;
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



// function to run the simulation
int main() {
    int N = 100; // Grid size
    int M = 10; // Number of simulations
    double stepSize = 0.1;
    

    // Output file for CSV results
    std::ofstream csvFile("results.csv");
    csvFile << "p, Average Steps, Std Dev Steps, Average Time Elapsed, Std Dev Time Elapsed, Average Fire Reached Bottom, Std Dev Fire Reached Bottom" << std::endl;



    std::vector<SimulationResult> results;
    int fireReachedBottom = 0;
    int steps = 0;


    

    
    for (double p = 0.0; p <= 1.0; p += stepSize) {

        for (int m=0; m<M ; m++){

            std::vector<std::vector<int>> grid(N, std::vector<int>(N, 0));
            initialiseGrid(grid, p, N);
            fireReachedBottom = 0;
            // vary wind direction from 0 to 3 to obtain different directions
            int windDirection = 0;
            SimulationResult result = simulateForestFire(grid, N, fireReachedBottom, steps, windDirection);
            #pragma omp critical
            results.push_back(result);
            
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
