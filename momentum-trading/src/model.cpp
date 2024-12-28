#include <iostream>
#include <vector>
#include <cmath>
#include <limits>
#include "utils.h"
#include "datasets.h"


int main() {
    ADAAI::CalibrationDataset Calibration = ADAAI::CalibrationDataset();
    ADAAI::ValidationDataset Validation = ADAAI::ValidationDataset();

    // Parameters
    std::vector<double> W_vals = { 5, 10, 15, 20, 30, 60 }; // Minutes

    std::vector<double> best_q;
    double max_sortino_ratio = -std::numeric_limits<double>::infinity();

    // Grid Search with OpenACC
    #pragma acc parallel loop collapse(5) copyin(CalibrationDataset[0:CalibrationDataset.size()])
    for (double alpha_1 = -1.0; alpha_1 < 10.0; alpha_1 += 10e-3) {
        for (double alpha_2 = -1.0; alpha_2 < 10.0; alpha_2 += 10e-3) {
            for (double alpha_3 = -1.0; alpha_3 < 10.0; alpha_3 += 10e-3) {
                for (double alpha_4 = -1.0; alpha_4 < 10.0; alpha_4 += 10e-3) {
                    for (double a = 0; a <= 10; a += 10e-3) {
                        for (double b = 10e-6; b <= 3; b += 10e-3) {
                            for (double w_ind = 0; w_ind < W_vals.size(); w_ind += 1){
                                // New params
                                std::vector q = {
                                    alpha_1,
                                    alpha_2,
                                    alpha_3,
                                    alpha_4,
                                    a,
                                    b,
                                    W_vals[w_ind],
                                };

                                // Sortino ratio computing
                                double avg_sortino_ratio = computeSortinoRatio(q, Calibration);

                                // Update best result
                                #pragma acc atomic update
                                if (avg_sortino_ratio > max_sortino_ratio) {
                                    max_sortino_ratio = avg_sortino_ratio;
                                    best_q = q;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    double sigma = 0.0; // Volatility measured over the whole validation set
    double T = Validation.size(); // Temporal length of the validation set

    double sortino_ratio = computeSortinoRatio(best_q, Validation);
    double pnl = computePnL(best_q, Validation);

    // Acceptance tests
    if (sortino_ratio > 1.0 &&
        (pnl > 2 * sigma * sqrt(T) || pnl > Validation.get_buy_and_hold_pnl() )) {
        std::cout << "Acceptance criteria passed!" << std::endl;
        }

    // Results
    std::cout << "Best parameters: ";
    for (const auto& param : best_q) {
        std::cout << param << " ";
    }
    std::cout << "\nMax Sortino Ratio: " << max_sortino_ratio << std::endl;

    return 0;

}
