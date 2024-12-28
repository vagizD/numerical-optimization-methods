#pragma once

namespace ADAAI {
    class Dataset {
    public:
        Dataset() {

        }
        size_t size();

        double get_buy_and_hold_pnl();
    };
    class CalibrationDataset : public Dataset {
    public:
        CalibrationDataset(): Dataset() {
        }

        size_t size();
        double get_buy_and_hold_pnl();
    };
    class ValidationDataset : public Dataset {
    public:
        ValidationDataset(): Dataset() {

        }

        size_t size();
        double get_buy_and_hold_pnl();
    };
}
