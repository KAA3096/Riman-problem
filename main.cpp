#include <iostream>
#include <cmath>
#include <cfloat>
#include <vector>
#include <fstream>
#include <algorithm>

const double PI = 3.14159265358979323846;

double f(double q) {
    return 0.5 * q * q;// нелинейный случай
    //return 1 * q;// линейный случай
}

double a(double q) {
    return q;
}

double minmod(double a, double b, double c) {
    if (a * b <= 0 || b * c <= 0) return 0.0;
    return std::copysign(1.0, a) * std::min({ std::abs(a), std::abs(b), std::abs(c) });
}

void solve(const bool isFirst, const double leftV, const double rightV, const double T, const double X, const int N, const double eps, const double CFL, const bool isSin) {

	std::vector<double> q(N, 0);
	std::vector<double> q_i_l(N + 1);
    std::vector<double> q_i_r(N + 1);
    std::vector<double> q_new(N);

    std::vector<double> w(N + 1, 0);
    std::vector<double> x(N + 1);
    std::vector<double> h(N);
    std::vector<double> S(N);

    std::vector<double> a_l(N + 1);
    std::vector<double> a_r(N + 1);

    std::vector<double> lambda_l(N + 1);
    std::vector<double> lambda_r(N + 1);

    std::vector<double> q_azr(N + 1);
    std::vector<double> q_ozr(N + 1);

    std::vector<double> q_i(N + 1);

    std::vector<double> F(N + 1);

    double tau;
    double t = 0;
    double dx = X / (N - 2);


    for (int i = 0; i <= N; i++)
        x[i] = (i - 1) * dx;
    for (int i = 0; i < N; i++)
        h[i] = dx;

    double target_x = X / 2;
    int break_index = 0;
    for (int i = 0; i < N - 1; ++i) {
        double x_center = 0.5 * (x[i] + x[i + 1]);
        if (x_center >= target_x) {
            break_index = i;
            break;
        }
    }

    const double qL = leftV;
    const double qR = rightV;
    for (int i = 1; i < N-1; ++i) {
        if (isSin) {
            double x_center = 0.5 * (x[i] + x[i + 1]);
            q[i] = std::sin(2.0 * PI * x_center / 1.0);
        }
        else {
			if (i < N / 2)
				q[i] = qL;
			else if (i > N / 2)
				q[i] = qR;
			else
				q[i] = (qL + qR) * 0.5;
        }
    }

    if (!isSin)
        q[N - 1] = qR;

    if (isFirst) {
        for (int i = 1; i < N - 1; i++)
            S[i] = 0;
    }
    else {
        for (int i = 0; i < N - 1; i++)
            S[i] = (q[i + 1] - q[i]) / h[i];
    }
    S[N - 1] = 0;

    while (t < T) {
        //Шаг 1
        {
            for (int i = 1; i < N; i++) {
                q_i_l[i] = q[i - 1] + 0.5 * h[i - 1] * S[i - 1];
                q_i_r[i] = q[i] - 0.5 * h[i] * S[i];
            }
        }

        //Шаг 2
        {
            for (int i = 1; i < N; i++) {
                a_l[i] = a(q_i_l[i]);
                a_r[i] = a(q_i_r[i]);
                double lambda_s = std::abs(q_i_l[i] - q_i_r[i]) > eps ? (f(q_i_r[i]) - f(q_i_l[i])) / (q_i_r[i] - q_i_l[i]) : a_l[i];

                lambda_l[i] = std::min((lambda_s), (a_l[i]));
                lambda_r[i] = std::max((lambda_s), (a_r[i]));
            }
        }

        //Шаг 3
        {
            double max_speed = 0.0;
            for (int i = 0; i < N; i++) {
                double abs_speed = std::abs(a(q[i])); // a(q) = q для f=0.5*q²
                if (abs_speed > max_speed) max_speed = abs_speed;
            }
            tau = CFL * dx / max_speed;
            if (t + tau > T) tau = T - t;
        }

        //Шаг 4
        {
            for (int i = 1; i < N; i++) {
                if (w[i] - lambda_l[i] <= eps) {
                    q_azr[i] = q_i_l[i];
                }
                else if (w[i] - lambda_r[i] >= eps ) {
                    q_azr[i] = q_i_r[i];
                }
                else {
                    q_azr[i] = 0;// 1.0 / a(w[i]);
                }
            }

            for (int i = 1; i < N; i++) {
                if (w[i] - lambda_l[i] <= eps) {
                    q_ozr[i] = S[i - 1] * (w[i] - a_l[i]);
                }
                else if (w[i] - lambda_r[i] >= eps) {
                    q_ozr[i] = S[i] * (w[i] - a_r[i]);
                }
                else {
                    q_ozr[i] = 0;
                }
            }
        }

        //Шаг 5
        {
            for (int i = 1; i < N; i++) {
                F[i] = q_azr[i] * w[i] - f(q_azr[i]) + 0.5 * q_ozr[i] * (w[i] - a(q_azr[i])) * tau;
            }

            for (int i = 1; i < N - 1; i++) {
                double h_new = x[i + 1] - x[i];
                q_new[i] = h[i] / h_new * q[i] + tau / h_new * (F[i + 1] - F[i]);
            }
            q_new[0] = q[0];
            q_new[N - 1] = q[N - 1];
        }

        //Шаг 6
        {
            for (int i = 1; i < N - 1; i++) {
                q_i[i] = q_azr[i] + q_ozr[i] * tau;
            }

            if (isFirst) {
                for (int i = 1; i < N - 1; i++)
                    S[i] = 0;
            }
            else {
                for (int i = 1; i < N - 1; i++) {
                    int prev = i - 1;// ( i == 0 ) ? N - 1 : i;
                    int next = i + 1;// ( i == N - 1 ) ? 0 : i + 1;
                    double S_raw = (q_i[i + 1] - q_i[i]) / h[i];
                    double dq_prev = q_i[i] - q_i[prev];
                    double dq_next = q_i[next] - q_i[i];
                    S[i] = minmod(dq_prev, h[i] * S_raw, dq_next) / h[i];
                }
            }


        }

        for (int i = 0; i < N; i++) {
            w[i] = (x[i + 1] - x[i]) / tau;
            h[i] = x[i + 1] - x[i];
            w[i] = 0;
            h[i] = dx;
        }
        q = q_new;

        t += tau;
    }

    //численное решение
    std::ofstream out(isFirst ? "results_G1.txt" : "results_G2.txt");
    out.precision(10);
    for (int i = 0; i < N; ++i) {
        double x_center = 0.5 * (x[i] + x[i + 1]);
        x_center = x_center < 0 ? 0 : x_center;
        x_center = x_center > X ? X : x_center;

        out << x_center << "\t" << q[i] << std::endl;
    }
    out.close();

    //точное решение
    std::ofstream outExact("resultsExact.txt");
    outExact.precision(10);
    for (int i = 0; i < N; ++i) {
        double x_center = 0.5 * (x[i] + x[i + 1]);
        x_center = x_center < 0 ? 0 : x_center;
        x_center = x_center > X ? X : x_center;

        double q_exact = 0;
        if (isSin) {
            q_exact = sin(2 * PI * (x_center - q[i] * t) / 1.0);
        }
        else {
            double x_exact_pos = X / 2 + 0.5 * t;
            q_exact = (x_center < x_exact_pos) ? qL : qR;
        }

        outExact << x_center << "\t" << q_exact << std::endl;
    }
    outExact.close();
}

int main() {
    //для синуса
    //const bool isFirst = false; //true - классическая схема, false - модифицированная схема
    const double leftV = 0.0;
    const double rightV = 1.0;

    const double T = 1.0;
    const double X = 1;
    const int N = 52; //на 2 больше чтобы впихнуть границы слева и справа
    const double eps = 0.000001;
    const double CFL = 0.5;
    const bool isSin = true; //true - синус, false - ступенька

    //для ступеньки
    //const bool isFirst = false; //true - классическая схема, false - модифицированная схема
   /* const double leftV = 0.0;
    const double rightV = 1.0;

    const double T = 1.0;
    const double X = 10;
    const int N = 78; //на 2 больше чтобы впихнуть границы слева и справа
    const double eps = 0.000001;
    const double CFL = 0.8;
    const bool isSin = false; //true - синус, false - ступенька*/

    solve(true, leftV, rightV, T, X, N, eps, CFL, isSin); //классическая схема 
    solve(false, leftV, rightV, T, X, N, eps, CFL, isSin); //модифицированная схема
}
