#include <iostream>
using namespace std;
/*
y1''(t) = f1(t, y1(t), y2(t), ..., yn(t), y1'(t), y2'(t), ..., yn'(t))
y2''(t) = f2(t, y1(t), y2(t), ..., yn(t), y1'(t), y2'(t), ..., yn'(t))
...
yn''(t) = fn(t, y1(t), y2(t), ..., yn(t), y1'(t), y2'(t), ..., yn'(t))

yi'(t) = vi(t)
vi'(t) = fi(t, y1(t), y2(t), ..., yn(t), v1(t), v2(t), ..., vn(t))
 */

void integrate(void (*f)(const double* y, double* k, double t, double* v), double t, double h, int N, double* y, double* v){

    double k1y[N], k2y[N], k3y[N], k4y[N], k5y[N], k6y[N], k7y[N], temp[N];
    double k1v[N], k2v[N], k3v[N], k4v[N], k5v[N], k6v[N], k7v[N];

    // Calculate k1v
    for (int i=0; i < N; i++) {
        k1y[i] = v[i];
    }
    f(y, k1v, t, k1y);

    // Calculate k2v

    for (int i = 0; i < N; i++) {
        temp[i] = y[i] + (1.0 / 5.0) * h * k1v[i];
        k2y[i] = v[i] +  (1.0 / 5.0) * h * k1v[i];
    }
    f(temp, k2v, t + (1.0 / 5.0) * h, k2y);

    // Calculate k3v

    for (int i = 0; i < N; i++) {
        temp[i] = y[i] + (3.0 / 40.0) * h * k1v[i] + (9.0 / 40.0) * h * k2v[i];
        k3y[i] = v[i] + (3.0 / 40.0) * h * k1v[i] + (9.0 / 40.0) * h * k2v[i];
    }
    f(temp, k3v, t + (3.0 / 10.0) * h, k3y);

    // Calculate k4v
    for (int i = 0; i < N; i++) {
        temp[i] = y[i] + (44.0 / 45.0) * h * k1v[i] - (56.0 / 15.0) * h * k2v[i] + (32.0 / 9.0) * h * k3v[i];
        k4y[i] = v[i] + (44.0 / 45.0) * h * k1v[i] - (56.0 / 15.0) * h * k2v[i] + (32.0 / 9.0) * h * k3v[i];
    }
    f(temp, k4v, t + (4.0 / 5.0) * h, k4y);

    // Calculate k5v
    for (int i = 0; i < N; i++) {
        temp[i] = y[i] + (19372.0 / 6561.0) * h * k1v[i] - (25360.0 / 2187.0) * h * k2v[i] +
                  (64448.0 / 6561.0) * h * k3v[i] - (212.0 / 729.0) * h * k4v[i];
        k5y[i] = v[i] + (19372.0 / 6561.0) * h * k1v[i] - (25360.0 / 2187.0) * h * k2v[i] +
                  (64448.0 / 6561.0) * h * k3v[i] - (212.0 / 729.0) * h * k4v[i];
    }
    f(temp, k5v, t + (8.0 / 9.0) * h, k5y);

    // Calculate k6v
    for (int i = 0; i < N; i++) {
        temp[i] = y[i] + (9017.0 / 3168.0) * h * k1v[i] - (355.0 / 33.0) * h * k2v[i] +
                  (46732.0 / 5247.0) * h * k3v[i] + (49.0 / 176.0) * h * k4v[i] - (5103.0 / 18656.0) * h * k5v[i];
        k6y[i] = v[i] + (9017.0 / 3168.0) * h * k1v[i] - (355.0 / 33.0) * h * k2v[i] +
                  (46732.0 / 5247.0) * h * k3v[i] + (49.0 / 176.0) * h * k4v[i] - (5103.0 / 18656.0) * h * k5v[i];
    }
    f(temp, k6v, t + h, k6y);

    // Calculate k7v
    for (int i = 0; i < N; i++) {
        temp[i] = y[i] + (35.0 / 384.0) * h * k1v[i] + (500.0 / 1113.0) * h * k3v[i] +
                  (125.0 / 192.0) * h * k4v[i] - (2187.0 / 6784.0) * h * k5v[i] + (11.0 / 84.0) * h * k6v[i];
        k7y[i] = v[i] + (35.0 / 384.0) * h * k1v[i] + (500.0 / 1113.0) * h * k3v[i] +
                  (125.0 / 192.0) * h * k4v[i] - (2187.0 / 6784.0) * h * k5v[i] + (11.0 / 84.0) * h * k6v[i];
    }
    f(temp, k7v, t + h, k7y);

    // Update state
    for (int i = 0; i < N; i++) {
        v[i] += (35.0 / 384.0) * k1v[i] + (500.0 / 1113.0) * k3v[i] + (125.0 / 192.0) * k4v[i] -
                (2187.0 / 6784.0) * k5v[i] + (11.0 / 84.0) * k6v[i];
        y[i] += (35.0 / 384.0) * k1y[i] + (500.0 / 1113.0) * k3y[i] + (125.0 / 192.0) * k4y[i] -
                (2187.0 / 6784.0) * k5y[i] + (11.0 / 84.0) * k6y[i];
    }

}

void derivative(const double* y, double* kv, double t, double* ky){
    kv[0] = y[1] + ky[0] * -0.5;
    kv[1] = -1.5 * y[0] + ky[1];
}

int main(){
    int N = 2;
    auto *y = new double [N];
    y[0] = 1.0; y[1] = 0.5;
    auto *v = new double [N];
    v[0] = 0.6; v[1] = -1.5;
    double t = 0.0, h = 0.01;
    for (int i = 0; i < 100; i++) {
        integrate(derivative, t, h, N, y, v);
        t += h;
        cout << y[0] << " " << y[1] << " " << v[0] << " " << v[1] <<endl;
    }

    return 0;
}
