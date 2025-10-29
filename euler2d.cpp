#include <iostream>
#include <iomanip>
#include <vector>
#include <array>
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>

enum SCHEME { CONSTANT = 0, MUSCL, WENO};
enum LIMITER { NONE = 0, MINMOD, VANLEER, VANALBADA};
enum FACE {WEST = 0, EAST, SOUTH, NORTH};

struct caseParameters {
  int numberOfPointsX;
  int numberOfPointsY;
  double gamma;
  double domainLengthX;
  double domainLengthY;
  double endTime;
  double CFL;
  double dx;
  double dy;
  double time;
  int timeStep;
};

int main() {

  /* ---------- Pre-processing ---------- */
    
    auto numericalScheme = SCHEME::MUSCL;
    auto limiter = LIMITER::VANLEER;
    caseParameters parameters;

    parameters.numberOfPointsX = 101;
    parameters.numberOfPointsY = 101;
    parameters.gamma = 1.4;
    parameters.domainLengthX = 1.0;
    parameters.domainLengthY = 1.0;
    parameters.endTime = 0.2;
    parameters.CFL = 0.1;
    
    parameters.dx = parameters.domainLengthX / (parameters.numberOfPointsX - 1);
    parameters.dy = parameters.domainLengthY / (parameters.numberOfPointsY - 1);
    parameters.time = 0.0;
    parameters.timeStep = 0;

    // 2D arrays: [i][j] indexing
    std::vector<std::vector<double>> x(parameters.numberOfPointsX, std::vector<double>(parameters.numberOfPointsY));
    std::vector<std::vector<double>> y(parameters.numberOfPointsX, std::vector<double>(parameters.numberOfPointsY));
    
    // U has 4 conservative variables in 2D: [rho, rho*u, rho*v, E]
    std::vector<std::vector<std::array<double, 4>>> U(parameters.numberOfPointsX, 
                                                       std::vector<std::array<double, 4>>(parameters.numberOfPointsY));
    
    // Ufaces: [i][j][face][variable] - face can be WEST, EAST, SOUTH, NORTH
    std::vector<std::vector<std::array<std::array<double, 4>, 4>>> Ufaces(parameters.numberOfPointsX,
                                                                            std::vector<std::array<std::array<double, 4>, 4>>(parameters.numberOfPointsY));
    
    // Ffaces for x-direction fluxes and y-direction fluxes
    std::vector<std::vector<std::array<std::array<double, 4>, 4>>> Ffaces(parameters.numberOfPointsX,
                                                                            std::vector<std::array<std::array<double, 4>, 4>>(parameters.numberOfPointsY));

    // Initialize mesh
    for (int i = 0; i < parameters.numberOfPointsX; i++) {
      for (int j = 0; j < parameters.numberOfPointsY; j++) {
        x[i][j] = i * parameters.dx;
        y[i][j] = j * parameters.dy;
      }
    }

    // Initial conditions - 2D Riemann problem (quadrant configuration)
    double rho = 0.0;
    double u = 0.0;
    double v = 0.0;
    double p = 0.0;

    for (int i = 0; i < parameters.numberOfPointsX; i++) {
      for (int j = 0; j < parameters.numberOfPointsY; j++) {
        // Simple 2D shock tube: high pressure in lower-left quadrant
        if (x[i][j] <= 0.5 && y[i][j] <= 0.5) {
          rho = 1.0;
          u = 0.0;
          v = 0.0;
          p = 1.0; 
        } else {
          rho = 0.125;
          u = 0.0;
          v = 0.0;
          p = 0.1; 
        }

        U[i][j][0] = rho;
        U[i][j][1] = rho * u;
        U[i][j][2] = rho * v;
        U[i][j][3] = p / (parameters.gamma - 1.0) + 0.5 * rho * (std::pow(u, 2) + std::pow(v, 2));
      }
    }

  /* ---------- Solving ---------- */
    while (parameters.time < parameters.endTime) {

      auto UOld = U;

      // calculate stable time step
      double speedMax = 0.0;
      for (int i = 0; i < parameters.numberOfPointsX; i++) {
        for (int j = 0; j < parameters.numberOfPointsY; j++) {
          auto rho = U[i][j][0];
          auto u = U[i][j][1] / rho;
          auto v = U[i][j][2] / rho;
          auto p = (parameters.gamma - 1.0) * (U[i][j][3] - 0.5 * rho * (std::pow(u, 2) + std::pow(v, 2)));
          
          // calculate wave speed for each cell
          double speedOfSound = std::sqrt(parameters.gamma * p / rho);
          double speedTotal = speedOfSound + std::sqrt(std::pow(u, 2) + std::pow(v, 2));
          if (speedTotal > speedMax)
            speedMax = speedTotal;
        }
      }
      double dt = (parameters.CFL * std::min(parameters.dx, parameters.dy)) / speedMax;

      // Solve equations - Reconstruction step

      if (numericalScheme == SCHEME::CONSTANT) {
        for (int i = 0; i < parameters.numberOfPointsX; i++) {
          for (int j = 0; j < parameters.numberOfPointsY; j++) {
            for (int variable = 0; variable < 4; ++variable) {
              Ufaces[i][j][FACE::WEST][variable] = U[i][j][variable];
              Ufaces[i][j][FACE::EAST][variable] = U[i][j][variable];
              Ufaces[i][j][FACE::SOUTH][variable] = U[i][j][variable];
              Ufaces[i][j][FACE::NORTH][variable] = U[i][j][variable];
            }
          }
        }
      } else if (numericalScheme == SCHEME::MUSCL) {
        // X-direction reconstruction
        for (int j = 0; j < parameters.numberOfPointsY; j++) {
          // use lower-order scheme near x-boundaries
          for (int variable = 0; variable < 4; ++variable) {
            Ufaces[0][j][FACE::WEST][variable] = U[0][j][variable];
            Ufaces[0][j][FACE::EAST][variable] = U[0][j][variable];

            Ufaces[parameters.numberOfPointsX - 1][j][FACE::WEST][variable] = U[parameters.numberOfPointsX - 1][j][variable];
            Ufaces[parameters.numberOfPointsX - 1][j][FACE::EAST][variable] = U[parameters.numberOfPointsX - 1][j][variable];
          }

          // use high-resolution MUSCL scheme on interior nodes / cells in x-direction
          for (int i = 1; i < parameters.numberOfPointsX - 1; i++) {
            for (int variable = 0; variable < 4; ++variable) {
              auto du_i_plus_half = U[i + 1][j][variable] - U[i][j][variable];
              auto du_i_minus_half = U[i][j][variable] - U[i - 1][j][variable];
              
              double rL = du_i_minus_half / (du_i_plus_half + 1e-8);
              double rR = du_i_plus_half / (du_i_minus_half + 1e-8);

              double psiL = 1.0;
              double psiR = 1.0;
              
              // apply limiter
              if (limiter == LIMITER::MINMOD) {
                psiL = std::max(0.0, std::min(1.0, rL));
                psiR = std::max(0.0, std::min(1.0, rR));
              } else if (limiter == LIMITER::VANLEER) {
                psiL = (rL + std::fabs(rL)) / (1.0 + std::fabs(rL));
                psiR = (rR + std::fabs(rR)) / (1.0 + std::fabs(rR));
              } else if (limiter == LIMITER::VANALBADA) {
                psiL = (rL * rL + rL) / (rL * rL + 1.0);
                psiR = (rR * rR + rR) / (rR * rR + 1.0);
              }

              Ufaces[i][j][FACE::WEST][variable] = U[i][j][variable]
                - 0.5 * psiL * du_i_plus_half;
              Ufaces[i][j][FACE::EAST][variable] = U[i][j][variable]
                + 0.5 * psiR * du_i_minus_half;
            }
          }
        }

        // Y-direction reconstruction
        for (int i = 0; i < parameters.numberOfPointsX; i++) {
          // use lower-order scheme near y-boundaries
          for (int variable = 0; variable < 4; ++variable) {
            Ufaces[i][0][FACE::SOUTH][variable] = U[i][0][variable];
            Ufaces[i][0][FACE::NORTH][variable] = U[i][0][variable];

            Ufaces[i][parameters.numberOfPointsY - 1][FACE::SOUTH][variable] = U[i][parameters.numberOfPointsY - 1][variable];
            Ufaces[i][parameters.numberOfPointsY - 1][FACE::NORTH][variable] = U[i][parameters.numberOfPointsY - 1][variable];
          }

          // use high-resolution MUSCL scheme on interior nodes / cells in y-direction
          for (int j = 1; j < parameters.numberOfPointsY - 1; j++) {
            for (int variable = 0; variable < 4; ++variable) {
              auto dv_j_plus_half = U[i][j + 1][variable] - U[i][j][variable];
              auto dv_j_minus_half = U[i][j][variable] - U[i][j - 1][variable];
              
              double rL = dv_j_minus_half / (dv_j_plus_half + 1e-8);
              double rR = dv_j_plus_half / (dv_j_minus_half + 1e-8);

              double psiL = 1.0;
              double psiR = 1.0;
              
              // apply limiter
              if (limiter == LIMITER::MINMOD) {
                psiL = std::max(0.0, std::min(1.0, rL));
                psiR = std::max(0.0, std::min(1.0, rR));
              } else if (limiter == LIMITER::VANLEER) {
                psiL = (rL + std::fabs(rL)) / (1.0 + std::fabs(rL));
                psiR = (rR + std::fabs(rR)) / (1.0 + std::fabs(rR));
              } else if (limiter == LIMITER::VANALBADA) {
                psiL = (rL * rL + rL) / (rL * rL + 1.0);
                psiR = (rR * rR + rR) / (rR * rR + 1.0);
              }

              Ufaces[i][j][FACE::SOUTH][variable] = U[i][j][variable]
                - 0.5 * psiL * dv_j_plus_half;
              Ufaces[i][j][FACE::NORTH][variable] = U[i][j][variable]
                + 0.5 * psiR * dv_j_minus_half;
            }
          }
        }
      }      

      // Compute fluxes at faces in X-direction
      std::array<double, 4> fluxL, fluxR;
      for (int i = 1; i < parameters.numberOfPointsX - 1; i++) {
        for (int j = 1; j < parameters.numberOfPointsY - 1; j++) {
          // WEST and EAST faces (x-direction fluxes)
          for (int face = FACE::WEST; face <= FACE::EAST; ++face) {
            int indexOffset = 0;
            if (face == FACE::WEST) indexOffset = 0;
            else if (face == FACE::EAST) indexOffset = 1;
          
            auto rhoL = Ufaces[i - 1 + indexOffset][j][FACE::EAST][0];
            auto uL = Ufaces[i - 1 + indexOffset][j][FACE::EAST][1] / rhoL;
            auto vL = Ufaces[i - 1 + indexOffset][j][FACE::EAST][2] / rhoL;
            auto EL = Ufaces[i - 1 + indexOffset][j][FACE::EAST][3];
            auto pL = (parameters.gamma - 1.0) * (EL - 0.5 * rhoL * (std::pow(uL, 2) + std::pow(vL, 2)));
            auto aL = std::sqrt(parameters.gamma * pL / rhoL);

            auto rhoR = Ufaces[i + indexOffset][j][FACE::WEST][0];
            auto uR = Ufaces[i + indexOffset][j][FACE::WEST][1] / rhoR;
            auto vR = Ufaces[i + indexOffset][j][FACE::WEST][2] / rhoR;
            auto ER = Ufaces[i + indexOffset][j][FACE::WEST][3];
            auto pR = (parameters.gamma - 1.0) * (ER - 0.5 * rhoR * (std::pow(uR, 2) + std::pow(vR, 2)));
            auto aR = std::sqrt(parameters.gamma * pR / rhoR);

            // X-direction flux
            fluxL[0] = rhoL * uL;
            fluxL[1] = pL + rhoL * std::pow(uL, 2);
            fluxL[2] = rhoL * uL * vL;
            fluxL[3] = uL * (EL + pL);

            fluxR[0] = rhoR * uR;
            fluxR[1] = pR + rhoR * std::pow(uR, 2);
            fluxR[2] = rhoR * uR * vR;
            fluxR[3] = uR * (ER + pR);

            // Rusanov Riemann solver
            auto speedMax = std::max(std::fabs(uL) + aL, std::fabs(uR) + aR);
            for (int variable = 0; variable < 4; ++variable) {
              const auto &qL = Ufaces[i - 1 + indexOffset][j][FACE::EAST][variable];
              const auto &qR = Ufaces[i + indexOffset][j][FACE::WEST][variable];
              const auto &fL = fluxL[variable];
              const auto &fR = fluxR[variable];
              Ffaces[i][j][face][variable] = 0.5 * (fL + fR) - 0.5 * speedMax * (qR - qL);
            }
          }

          // SOUTH and NORTH faces (y-direction fluxes)
          for (int face = FACE::SOUTH; face <= FACE::NORTH; ++face) {
            int indexOffset = 0;
            if (face == FACE::SOUTH) indexOffset = 0;
            else if (face == FACE::NORTH) indexOffset = 1;
          
            auto rhoL = Ufaces[i][j - 1 + indexOffset][FACE::NORTH][0];
            auto uL = Ufaces[i][j - 1 + indexOffset][FACE::NORTH][1] / rhoL;
            auto vL = Ufaces[i][j - 1 + indexOffset][FACE::NORTH][2] / rhoL;
            auto EL = Ufaces[i][j - 1 + indexOffset][FACE::NORTH][3];
            auto pL = (parameters.gamma - 1.0) * (EL - 0.5 * rhoL * (std::pow(uL, 2) + std::pow(vL, 2)));
            auto aL = std::sqrt(parameters.gamma * pL / rhoL);

            auto rhoR = Ufaces[i][j + indexOffset][FACE::SOUTH][0];
            auto uR = Ufaces[i][j + indexOffset][FACE::SOUTH][1] / rhoR;
            auto vR = Ufaces[i][j + indexOffset][FACE::SOUTH][2] / rhoR;
            auto ER = Ufaces[i][j + indexOffset][FACE::SOUTH][3];
            auto pR = (parameters.gamma - 1.0) * (ER - 0.5 * rhoR * (std::pow(uR, 2) + std::pow(vR, 2)));
            auto aR = std::sqrt(parameters.gamma * pR / rhoR);

            // Y-direction flux
            fluxL[0] = rhoL * vL;
            fluxL[1] = rhoL * uL * vL;
            fluxL[2] = pL + rhoL * std::pow(vL, 2);
            fluxL[3] = vL * (EL + pL);

            fluxR[0] = rhoR * vR;
            fluxR[1] = rhoR * uR * vR;
            fluxR[2] = pR + rhoR * std::pow(vR, 2);
            fluxR[3] = vR * (ER + pR);

            // Rusanov Riemann solver
            auto speedMax = std::max(std::fabs(vL) + aL, std::fabs(vR) + aR);
            for (int variable = 0; variable < 4; ++variable) {
              const auto &qL = Ufaces[i][j - 1 + indexOffset][FACE::NORTH][variable];
              const auto &qR = Ufaces[i][j + indexOffset][FACE::SOUTH][variable];
              const auto &fL = fluxL[variable];
              const auto &fR = fluxR[variable];
              Ffaces[i][j][face][variable] = 0.5 * (fL + fR) - 0.5 * speedMax * (qR - qL);
            }
          }
        }
      }

      // Update solution
      for (int i = 1; i < parameters.numberOfPointsX - 1; i++) {
        for (int j = 1; j < parameters.numberOfPointsY - 1; j++) {
          for (int k = 0; k < 4; k++) {
            const auto &dFx = Ffaces[i][j][FACE::EAST][k] - Ffaces[i][j][FACE::WEST][k];
            const auto &dFy = Ffaces[i][j][FACE::NORTH][k] - Ffaces[i][j][FACE::SOUTH][k];
            U[i][j][k] = UOld[i][j][k] - (dt / parameters.dx) * dFx - (dt / parameters.dy) * dFy;
          }
        }
      }

      // Update boundary conditions (simple extrapolation)
      for (int j = 0; j < parameters.numberOfPointsY; j++) {
        for (int k = 0; k < 4; k++) {
          U[0][j][k] = U[1][j][k];
          U[parameters.numberOfPointsX - 1][j][k] = U[parameters.numberOfPointsX - 2][j][k];
        }
      }
      for (int i = 0; i < parameters.numberOfPointsX; i++) {
        for (int k = 0; k < 4; k++) {
          U[i][0][k] = U[i][1][k];
          U[i][parameters.numberOfPointsY - 1][k] = U[i][parameters.numberOfPointsY - 2][k];
        }
      }

      // Output to file (VTK format for visualization)
      if (parameters.timeStep % 10 == 0) {
        std::ofstream outputFile;

        std::ostringstream timeStepTemp;
        timeStepTemp << std::setfill('0') << std::setw(6);
        timeStepTemp << parameters.timeStep;
        auto timeStep = timeStepTemp.str();

        outputFile.open("solution2d_" + timeStep + ".vtk");
        
        // VTK header
        outputFile << "# vtk DataFile Version 3.0\n";
        outputFile << "2D Euler Solution\n";
        outputFile << "ASCII\n";
        outputFile << "DATASET STRUCTURED_POINTS\n";
        outputFile << "DIMENSIONS " << parameters.numberOfPointsX << " " << parameters.numberOfPointsY << " 1\n";
        outputFile << "ORIGIN 0 0 0\n";
        outputFile << "SPACING " << parameters.dx << " " << parameters.dy << " 1\n";
        outputFile << "POINT_DATA " << parameters.numberOfPointsX * parameters.numberOfPointsY << "\n";
        
        // Density
        outputFile << "SCALARS density double\n";
        outputFile << "LOOKUP_TABLE default\n";
        for (int j = 0; j < parameters.numberOfPointsY; j++) {
          for (int i = 0; i < parameters.numberOfPointsX; i++) {
            outputFile << U[i][j][0] << "\n";
          }
        }
        
        // Pressure
        outputFile << "SCALARS pressure double\n";
        outputFile << "LOOKUP_TABLE default\n";
        for (int j = 0; j < parameters.numberOfPointsY; j++) {
          for (int i = 0; i < parameters.numberOfPointsX; i++) {
            auto rho = U[i][j][0];
            auto u = U[i][j][1] / rho;
            auto v = U[i][j][2] / rho;
            auto p = (parameters.gamma - 1.0) * (U[i][j][3] - 0.5 * rho * (std::pow(u, 2) + std::pow(v, 2)));
            outputFile << p << "\n";
          }
        }
        
        // Velocity
        outputFile << "VECTORS velocity double\n";
        for (int j = 0; j < parameters.numberOfPointsY; j++) {
          for (int i = 0; i < parameters.numberOfPointsX; i++) {
            auto rho = U[i][j][0];
            auto u = U[i][j][1] / rho;
            auto v = U[i][j][2] / rho;
            outputFile << u << " " << v << " 0\n";
          }
        }
        
        outputFile.close();
      }

      // output current time step information to screen
      std::cout << "Current time: " << std::scientific << std::setw(10) << std::setprecision(3) << parameters.time;
      std::cout << ", End time: " << std::scientific << std::setw(10) << std::setprecision(3) << parameters.endTime;
      std::cout << ", Current time step: " << std::fixed << std::setw(7) << parameters.timeStep;
      std::cout << "\r" << std::flush;

      parameters.time += dt;
      parameters.timeStep++;
    }
  
  std::cout << "\nSimulation finished !! :) :)" << std::endl;

  return 0;
}
