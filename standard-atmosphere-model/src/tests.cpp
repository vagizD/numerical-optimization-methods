#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <cmath>
#include <iostream>
#include <string>
#include "constants.h"
#include "integrator.h"
#include "observer.h"
#include "sam.h"
#include "timestepper.h"

int func(double t, const double y[], double f[], void *params) {
    (void)(t);
    const Pressure m_pressure(101325, 288.15);
    const Density m_density(m_pressure);
    const AirSpeed m_airSpeed(m_density);
    DragCoef m_cd;
    AerodynamicDragForce m_q(m_cd, m_density);
    double s = M_PI * 0.216 * 0.216 / 4;
    double velocity = std::sqrt(y[1] * y[1] + y[3] * y[3]);
    double mach = m_airSpeed(y[2]);
    double qdf = m_q(mach, y[2], velocity, s);
    f[0] = y[1];
    f[1] = -(qdf * std::abs(y[1])) / (106.0 * velocity);
    f[2] = y[3];
    f[3] = -(qdf * std::abs(y[3])) / (106.0 * velocity) - g;
    return GSL_SUCCESS;
}

int jac(double t, const double y[], double *dfdy, double dfdt[], void *params) {
    (void)(t);
    const Pressure m_pressure(101325, 288.15);
    const Density m_density(m_pressure);
    const AirSpeed m_airSpeed(m_density);
    DragCoef m_cd;
    AerodynamicDragForce m_q(m_cd, m_density);
    double s = M_PI * 0.216 * 0.216 / 4;
    double velocity = std::sqrt(y[1] * y[1] + y[3] * y[3]);
    double mach = m_airSpeed(y[2]);
    double qdf = m_q(mach, y[2], velocity, s);
    gsl_matrix_view dfdy_mat = gsl_matrix_view_array(dfdy, 4, 4);
    gsl_matrix *m = &dfdy_mat.matrix;
    gsl_matrix_set(m, 0, 0, 0.0);
    gsl_matrix_set(m, 0, 1, 1.0);
    gsl_matrix_set(m, 0, 2, 0.0);
    gsl_matrix_set(m, 0, 3, 0.0);
    gsl_matrix_set(m, 1, 0, 0.0);
    gsl_matrix_set(m, 1, 1, -(qdf) / (106.0 * velocity));
    gsl_matrix_set(m, 1, 2, 0.0);
    gsl_matrix_set(m, 1, 3, 0.0);
    gsl_matrix_set(m, 2, 0, 0.0);
    gsl_matrix_set(m, 2, 1, 0.0);
    gsl_matrix_set(m, 2, 2, 0.0);
    gsl_matrix_set(m, 2, 3, 1.0);
    gsl_matrix_set(m, 3, 0, 0.0);
    gsl_matrix_set(m, 3, 1, 0.0);
    gsl_matrix_set(m, 3, 2, 0.0);
    gsl_matrix_set(m, 3, 3, -(qdf) / (106.0 * velocity) - g);
    dfdt[0] = 0.0;
    dfdt[1] = 0.0;
    dfdt[2] = 0.0;
    dfdt[3] = 0.0;
    return GSL_SUCCESS;
}

void save_projectile_gsl(std::string &filename, const std::vector<Point> &a_projectile) {
    std::ofstream fhandle;
    fhandle.open(filename);
    fhandle << "t,x,y\n";
    for (auto &p : a_projectile) {
        fhandle << p.get_t() << ", " << p.get_x() << ", " << p.get_y() << '\n';
    }
    fhandle.close();
}

int main() {
    double p0 = 101325;
    double t0 = 288.15;
    double mass = 106.0;
    double diameter = 0.216;

    double alpha = 5;  // degrees

    double max_range = -1;
    double max_range_alpha = -1;
    double max_height = -1;
    double max_height_alpha = -1;
    double max_range_gsl = -1;
    double max_range_alpha_gsl = -1;
    double max_height_gsl = -1;
    double max_height_alpha_gsl = -1;

    for (int counter = 25; counter <= 60; counter++) {
        std::cout << "Checking the alpha: " << alpha << std::endl;
        double v0 = 1640;  // meter/second
        double x0 = 0;
        double y0 = 0;

        std::string path = "../projectiles/trajectory/trajectory_optimal_range.csv";

        auto rhs = RHS_Projectile(p0, t0, mass, diameter);
        auto stepper = TimeStepper_RKF45<RHS_Projectile>(&rhs);
        auto observer = SimpleObserver(x0, y0);
        auto integrator =
            ODE_Integrator<TimeStepper_RKF45<RHS_Projectile>, SimpleObserver>(
                &stepper, &observer, path, false, false
            );

        double t_start = 0;
        double h_start = 1e-3;
        double h_gsl = h_start / 100.0;
        std::array<double, 4> u = {
            x0, v0 * std::cos(alpha * M_PI / 180), y0, v0 * std::sin(alpha * M_PI / 180)};
        std::array<double, 4> result = {};

        integrator(t_start, h_start, u, result);
        // std::cout << result[0] << std::endl;
        if (result[0] > max_range) {
            ;
            max_range = result[0];
            max_range_alpha = alpha;
        }
        if (observer.m_y_max > max_height) {
            max_height = observer.m_y_max;
            max_height_alpha = alpha;
        }

        // std::cout << "===== GSL impelementation =====\n";

        const gsl_odeiv2_step_type *T = gsl_odeiv2_step_rkf45;

        gsl_odeiv2_step *s = gsl_odeiv2_step_alloc(T, 4);
        gsl_odeiv2_control *c = gsl_odeiv2_control_y_new(1e-6, 0.0);
        gsl_odeiv2_evolve *e = gsl_odeiv2_evolve_alloc(4);

        double mu = 1;
        gsl_odeiv2_system sys = {func, jac, 4, &mu};

        double t = 0.0, t1 = 200.0;
        double y[4] = {
            x0, v0 * std::cos(alpha * M_PI / 180), y0, v0 * std::sin(alpha * M_PI / 180)};
        double prev_range = -1;
        while (y[2] >= 0) {
            if (max_height_gsl < y[2]) {
                max_height_gsl = y[2];
                max_height_alpha_gsl = alpha;
            }
            int status = gsl_odeiv2_evolve_apply(e, c, s, &sys, &t, t1, &h_gsl, y);

            if (status != GSL_SUCCESS) {
                std::cout << "error , return value=" << status << "\n";
                break;
            }
            if (y[2] < 0) {
                break;
            }
            prev_range = y[0];
            //        std::cout << "GSL-State: t = " << t << ", x = " << y[0]
            //                  << ", y = " << y[2] << std::endl;
        }
        if (prev_range > max_range_gsl) {
            max_range_gsl = prev_range;
            max_range_alpha_gsl = alpha;
        }
        gsl_odeiv2_evolve_free(e);
        gsl_odeiv2_control_free(c);
        gsl_odeiv2_step_free(s);
        alpha++;
    }

    std::cout << "Max range: " << max_range << "\nOptimal angle: " << max_range_alpha
              << "\n";
    std::cout << "Max height: " << max_height << "\nOptimal angle: " << max_height_alpha
              << "\n";
    std::cout << "=============================================================\n";
    std::cout << "Max range with GSL: " << max_range_gsl
              << "\nOptimal angle with GSL: " << max_range_alpha_gsl << "\n";
    std::cout << "Max height with GSL: " << max_height_gsl
              << "\nOptimal angle with GSL: " << max_height_alpha_gsl << "\n";

    {
        alpha = max_range_alpha;
        double v0 = 1640;  // meter/second
        double x0 = 0;
        double y0 = 0;

        std::string path = "../projectiles/trajectory/trajectory_optimal_range.csv";

        auto rhs = RHS_Projectile(p0, t0, mass, diameter);
        auto stepper = TimeStepper_RKF45<RHS_Projectile>(&rhs);
        auto observer = SimpleObserver(x0, y0);
        auto integrator =
            ODE_Integrator<TimeStepper_RKF45<RHS_Projectile>, SimpleObserver>(
                &stepper, &observer, path, false, true
            );

        double t_start = 0;
        double h_start = 1e-3;
        double h_gsl = h_start / 100.0;
        std::array<double, 4> u = {
            x0, v0 * std::cos(alpha * M_PI / 180), y0, v0 * std::sin(alpha * M_PI / 180)};
        std::array<double, 4> result = {};

        integrator(t_start, h_start, u, result);

        alpha = max_range_alpha_gsl;

        path = "../projectiles/trajectory/trajectory_optimal_range_gsl.csv";
        std::vector<Point> trajectory_gsl;

        const gsl_odeiv2_step_type *T = gsl_odeiv2_step_rkf45;

        gsl_odeiv2_step *s = gsl_odeiv2_step_alloc(T, 4);
        gsl_odeiv2_control *c = gsl_odeiv2_control_y_new(1e-6, 0.0);
        gsl_odeiv2_evolve *e = gsl_odeiv2_evolve_alloc(4);

        double mu = 1;
        gsl_odeiv2_system sys = {func, jac, 4, &mu};

        double t = 0.0, t1 = 200.0;
        double y[4] = {
            x0, v0 * std::cos(alpha * M_PI / 180), y0, v0 * std::sin(alpha * M_PI / 180)};
        while (y[2] >= 0) {
            trajectory_gsl.emplace_back(y[0], y[2], t);
            int status = gsl_odeiv2_evolve_apply(e, c, s, &sys, &t, t1, &h_gsl, y);

            if (status != GSL_SUCCESS) {
                std::cout << "error , return value=" << status << "\n";
                break;
            }
            if (y[2] < 0) {
                break;
            }
        }
        save_projectile_gsl(path, trajectory_gsl);
        gsl_odeiv2_evolve_free(e);
        gsl_odeiv2_control_free(c);
        gsl_odeiv2_step_free(s);
    }

    {
        alpha = max_height_alpha;
        double v0 = 1640;  // meter/second
        double x0 = 0;
        double y0 = 0;

        std::string path = "../projectiles/trajectory/trajectory_optimal_height.csv";

        auto rhs = RHS_Projectile(p0, t0, mass, diameter);
        auto stepper = TimeStepper_RKF45<RHS_Projectile>(&rhs);
        auto observer = SimpleObserver(x0, y0);
        auto integrator =
            ODE_Integrator<TimeStepper_RKF45<RHS_Projectile>, SimpleObserver>(
                &stepper, &observer, path, false, true
            );

        double t_start = 0;
        double h_start = 1e-3;
        double h_gsl = h_start / 100.0;
        std::array<double, 4> u = {
            x0, v0 * std::cos(alpha * M_PI / 180), y0, v0 * std::sin(alpha * M_PI / 180)};
        std::array<double, 4> result = {};

        integrator(t_start, h_start, u, result);

        alpha = max_height_alpha_gsl;

        path = "../projectiles/trajectory/trajectory_optimal_height_gsl.csv";
        std::vector<Point> trajectory_gsl;

        const gsl_odeiv2_step_type *T = gsl_odeiv2_step_rkf45;

        gsl_odeiv2_step *s = gsl_odeiv2_step_alloc(T, 4);
        gsl_odeiv2_control *c = gsl_odeiv2_control_y_new(1e-6, 0.0);
        gsl_odeiv2_evolve *e = gsl_odeiv2_evolve_alloc(4);

        double mu = 1;
        gsl_odeiv2_system sys = {func, jac, 4, &mu};

        double t = 0.0, t1 = 200.0;
        double y[4] = {
            x0, v0 * std::cos(alpha * M_PI / 180), y0, v0 * std::sin(alpha * M_PI / 180)};
        while (y[2] >= 0) {
            trajectory_gsl.emplace_back(y[0], y[2], t);
            int status = gsl_odeiv2_evolve_apply(e, c, s, &sys, &t, t1, &h_gsl, y);

            if (status != GSL_SUCCESS) {
                std::cout << "error , return value=" << status << "\n";
                break;
            }
            if (y[2] < 0) {
                break;
            }
        }
        save_projectile_gsl(path, trajectory_gsl);
        gsl_odeiv2_evolve_free(e);
        gsl_odeiv2_control_free(c);
        gsl_odeiv2_step_free(s);
    }
    return 0;
}
