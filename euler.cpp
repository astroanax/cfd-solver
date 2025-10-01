#include <iostream>
#include <iomanip>
#include <vector>
#include <array>
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>

enum SCHEME { CONSTANT = 0, MUSCL};
enum LIMITER { NONE = 0, MINMOD, VANLEER};
enum FACE {WEST = 0, EAST};

struct caseParameters {
  int numberOfPoints;
  double gamma;
  double domainLength;
  double endTime;
  double CFL;
  double dx;
  double time;
  int timeStep;
};

int main() {

  /* ---------- Pre-processing ---------- */
    
    auto numericalScheme = SCHEME::MUSCL;
    auto limiter = LIMITER::VANLEER;
    caseParameters parameters;

    parameters.numberOfPoints = 101;
    parameters.gamma = 1.4;
    parameters.domainLength = 1.0;
    parameters.endTime = 0.2;
    parameters.CFL = 0.1;
    
    parameters.dx = parameters.domainLength / (parameters.numberOfPoints - 1);
    parameters.time = 0.0;
    parameters.timeStep = 0;

    std::vector<double> x(parameters.numberOfPoints);
    std::vector<std::array<double, 3>> U(parameters.numberOfPoints);
    std::vector<std::array<std::array<double, 3>, 2>> Ufaces(parameters.numberOfPoints);
    std::vector<std::array<std::array<double, 3>, 2>> Ffaces(parameters.numberOfPoints);

    for (int i = 0; i < parameters.numberOfPoints; i++) {
      x[i] = i * parameters.dx;
    }

    double rho = 0.0;
    double u = 0.0;
    double p = 0.0;

    for (int i = 0; i < parameters.numberOfPoints; i++) {
      if (x[i] <= 0.5) {
        rho = 1.0;
        u = 0.0;
        p = 1.0; 
      } else {
        rho = 0.125;
        u = 0.0;
        p = 0.1; 
      }

      U[i][0] = rho;
      U[i][1] = rho * u;
      U[i][2] = p / (parameters.gamma - 1.0) + 0.5 * rho * std::pow(u, 2);
    }

  /* ---------- Solving ---------- */
    while (parameters.time < parameters.endTime) {

      auto UOld = U;

      // calculate stable time step
      double speedMax = 0.0;
      for (int i = 0; i < parameters.numberOfPoints; i++) {
        auto rho = U[i][0];
        auto u = U[i][1] / rho;
        auto p = (parameters.gamma - 1.0) * (U[i][2] - 0.5 * rho * std::pow(u, 2));
        
        // calculate wave speed for each cell
        double speedOfSound = std::sqrt(parameters.gamma * p / rho);
        if (speedOfSound + std::fabs(u) > speedMax)
          speedMax = speedOfSound + std::fabs(u);
      }
      double dt = (parameters.CFL * parameters.dx) / speedMax;

      // Solve equations

      if (numericalScheme == SCHEME::CONSTANT) {
        for (int i = 0; i < parameters.numberOfPoints; i++) {
          for (int variable = 0; variable < 3; ++variable) {
            Ufaces[i][FACE::WEST][variable] = U[i][variable];
            Ufaces[i][FACE::EAST][variable] = U[i][variable];
          }
        }
      } else if (numericalScheme == SCHEME::MUSCL) {
        // use lower-order scheme near boundaries
        for (int variable = 0; variable < 3; ++variable) {
          Ufaces[0][FACE::WEST][variable] = U[0][variable];
          Ufaces[0][FACE::EAST][variable] = U[0][variable];

          Ufaces[parameters.numberOfPoints - 1][FACE::WEST][variable] = U[parameters.numberOfPoints - 1][variable];
          Ufaces[parameters.numberOfPoints - 1][FACE::EAST][variable] = U[parameters.numberOfPoints - 1][variable];
        }

        // use high-resolution MUSCL scheme on interior nodes / cells
        for (int i = 1; i < parameters.numberOfPoints - 1; i++) {
          for (int variable = 0; variable < 3; ++variable) {
            auto du_i_plus_half = U[i + 1][variable] - U[i][variable];
            auto du_i_minus_half = U[i][variable] - U[i - 1][variable];
            
            double rL = du_i_minus_half / (du_i_plus_half + 1e-8);
            double rR = du_i_plus_half / (du_i_minus_half + 1e-8);

            double psiL = 1.0;
            double psiR = 1.0;
            
            // apply limiter to make scheme TVD (total variation diminishing)
            if (limiter == LIMITER::MINMOD) {
              psiL = std::max(0.0, std::min(1.0, rL));
              psiR = std::max(0.0, std::min(1.0, rR));
            } else if (limiter == LIMITER::VANLEER) {
              psiL = (rL + std::fabs(rL)) / (1.0 + std::fabs(rL));
              psiR = (rR + std::fabs(rR)) / (1.0 + std::fabs(rR));
            }

            Ufaces[i][FACE::WEST][variable] = U[i][variable]
              - 0.5 * psiL * du_i_plus_half;
            Ufaces[i][FACE::EAST][variable] = U[i][variable]
              + 0.5 * psiR * du_i_minus_half;
          }
        }
      }      

      // compute fluxes at faces
      std::array<double, 3> fluxL, fluxR;
      for (int i=1; i < parameters.numberOfPoints - 1; i++) {
        for (int face = FACE::WEST; face <= FACE::EAST; ++face) {
          int indexOffset = 0;
          if (face == FACE::WEST) indexOffset = 0;
          else if (face == FACE::EAST) indexOffset = 1;
        
          auto rhoL = Ufaces[i - 1 + indexOffset][FACE::EAST][0];
          auto uL = Ufaces[i - 1 + indexOffset][FACE::EAST][1] / rhoL;
          auto EL = Ufaces[i - 1 + indexOffset][FACE::EAST][2];
          auto pL = (parameters.gamma - 1.0) * (EL - 0.5 * rhoL * std::pow(uL, 2));
          auto aL = std::sqrt(parameters.gamma * pL / rhoL);

          auto rhoR = Ufaces[i + indexOffset][FACE::WEST][0];
          auto uR = Ufaces[i + indexOffset][FACE::WEST][1] / rhoR;
          auto ER = Ufaces[i + indexOffset][FACE::WEST][2];
          auto pR = (parameters.gamma - 1.0) * (ER - 0.5 * rhoR * std::pow(uR, 2));
          auto aR = std::sqrt(parameters.gamma * pR / rhoR);

          fluxL[0] = rhoL * uL;
          fluxL[1] = pL + rhoL * std::pow(uL, 2);
          fluxL[2] = uL * (EL + pL);

          fluxR[0] = rhoR * uR;
          fluxR[1] = pR + rhoR * std::pow(uR, 2);
          fluxR[2] = uR * (ER + pR);

          // Rusanov Riemann solver
          auto speedMax = std::max(std::fabs(uL) + aL, std::fabs(uR) + aR);
          for (int variable = 0; variable < 3; ++variable) {
            const auto &qL = Ufaces[i - 1 + indexOffset][FACE::EAST][variable];
            const auto &qR = Ufaces[i + indexOffset][FACE::WEST][variable];
            const auto &fL = fluxL[variable];
            const auto &fR = fluxR[variable];
            Ffaces[i][face][variable] = 0.5 * (fL + fR) - speedMax * (qR - qL);
          }
        }
      }
      for (int i=1; i < parameters.numberOfPoints - 1; i++)
        for (int j=0; j<3; j++) {
          const auto &dF = Ffaces[i][FACE::EAST][j] - Ffaces[i][FACE::WEST][j];
          U[i][j] = UOld[i][j] - (dt / parameters.dx) * dF;
        }

      // Update boundary conditions
      auto rhoL = 1.0;
      auto uL = 0.0;

      auto rhoR = 0.125;
      auto uR = 0.0;

      U[0][0] = U[1][0];
      U[0][1] = rhoL * uL;
      U[0][2] = U[1][2];
      
      U[parameters.numberOfPoints - 1][0] = U[parameters.numberOfPoints - 2][0];
      U[parameters.numberOfPoints - 1][1] = rhoR * uR;
      U[parameters.numberOfPoints - 1][2] = U[parameters.numberOfPoints - 2][2];

      std::ofstream outputFile;

      // convert time step and points into 6 digits string with leading zeros
      std::ostringstream timeStepTemp, pointsTemp;

      timeStepTemp << std::setfill('0') << std::setw(6);
      timeStepTemp << parameters.timeStep;
      auto timeStep = timeStepTemp.str();

      pointsTemp << std::setfill('0') << std::setw(6);
      pointsTemp << parameters.numberOfPoints;
      auto points = pointsTemp.str();

      outputFile.open("solution_" + points + "_" + timeStep + ".csv");
      outputFile << "x,rho,u,p" << std::endl;
      for (int i = 0; i < parameters.numberOfPoints; i++) {
        auto rho = U[i][0];
        auto u = U[i][1] / rho;
        auto p = (parameters.gamma - 1.0) * (U[i][2] - 0.5 * rho * std::pow(u, 2));
        outputFile << x[i] << "," << rho << "," << u << "," << p << std::endl;
      }
      outputFile.close();

      // output current time step information to screen
      std::cout << "Current time: " << std::scientific << std::setw(10) << std::setprecision(3) << parameters.time;
      std::cout << ", End time: " << std::scientific << std::setw(10) << std::setprecision(3) << parameters.endTime;
      std::cout << ", Current time step: " << std::fixed << std::setw(7) << parameters.timeStep;
      std::cout << "\r";

      parameters.time += dt;
      parameters.timeStep++;
    }
  
  std::cout << "\nSimulation finished !! :) :)" << std::endl;

  return 0;
}
