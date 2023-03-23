#include <iostream>
#include <cmath>
#include <vector>

using namespace std;

template <int N>
class DormandPrince {
public:
    typedef vector<double> state_type;

    typedef void (*deriv_func)(const state_type &state, state_type &derivative, double t);

    void integrate(deriv_func f, state_type &state, double t, double h) {
        state_type k1(N), k2(N), k3(N), k4(N), k5(N), k6(N), k7(N), temp(N);

        // Calculate k1
        f(state, k1, t);

        // Calculate k2
        for (int i = 0; i < N; i++) {
            temp[i] = state[i] + (1.0 / 5.0) * h * k1[i];
        }
        f(temp, k2, t + (1.0 / 5.0) * h);

        // Calculate k3
        for (int i = 0; i < N; i++) {
            temp[i] = state[i] + (3.0 / 40.0) * h * k1[i] + (9.0 / 40.0) * h * k2[i];
        }
        f(temp, k3, t + (3.0 / 10.0) * h);

        // Calculate k4
        for (int i = 0; i < N; i++) {
            temp[i] = state[i] + (44.0 / 45.0) * h * k1[i] - (56.0 / 15.0) * h * k2[i] + (32.0 / 9.0) * h * k3[i];
        }
        f(temp, k4, t + (4.0 / 5.0) * h);

        // Calculate k5
        for (int i = 0; i < N; i++) {
            temp[i] = state[i] + (19372.0 / 6561.0) * h * k1[i] - (25360.0 / 2187.0) * h * k2[i] +
                      (64448.0 / 6561.0) * h * k3[i] - (212.0 / 729.0) * h * k4[i];
        }
        f(temp, k5, t + (8.0 / 9.0) * h);

        // Calculate k6
        for (int i = 0; i < N; i++) {
            temp[i] = state[i] + (9017.0 / 3168.0) * h * k1[i] - (355.0 / 33.0) * h * k2[i] +
                      (46732.0 / 5247.0) * h * k3[i] + (49.0 / 176.0) * h * k4[i] - (5103.0 / 18656.0) * h * k5[i];
        }
        f(temp, k6, t + h);

        // Calculate k7
        for (int i = 0; i < N; i++) {
            temp[i] = state[i] + (35.0 / 384.0) * h * k1[i] + (500.0 / 1113.0) * h * k3[i] +
                      (125.0 / 192.0) * h * k4[i] - (2187.0 / 6784.0) * h * k5[i] + (11.0 / 84.0) * h * k6[i];
        }
        f(temp, k7, t + h);

        // Update state
        for (int i = 0; i < N; i++) {
            state[i] += (5179.0 / 57600.0) * k1[i] + (7571.0 / 16695.0) * k3[i] + (393.0 / 640.0) * k4[i] -
                        (92097.0 / 339200.0) * k5[i] + (187.0 / 2100.0) * k6[i] + (1.0 / 40.0) * k7[i];
        }
    }
};

void derivative(const DormandPrince<2>::state_type& state, DormandPrince<2>::state_type& derivative, double t)
{
    derivative[0] = state[1];
    derivative[1] = -1.0*state[0];
}

int main() {
    DormandPrince<2> dp;
    DormandPrince<2>::state_type state = {1.0, 0.0};
    double t = 0.0, h = 0.01;

    for (int i = 0; i < 100; i++) {
        dp.integrate(derivative, state, t, h);
        t += h;

        cout << state[0] << " " << state[1] << endl;
    }

    return 0;
}
